%==============================================================================
% This code is part of the Finite Element Method app for the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR/FAIRFEM 
%==============================================================================
%
% function [Sc,dS,d2S] = hyperElasticFEM(uc,yRef,Mesh,varargin)
%
% hyperelastic regularization based on tetrahedral finite element
% discretization with linear basis functions
%
% S(u)=\int a1*||\nabla u||^2 + a2*phi(cof\nabla(x+u)) + a3*psi(det\nabla(x+u))
%
% Input:
%   uc     - coefficients for displacement
%   yRef   - coefficients of reference transformation
%   Mesh   - description of mesh, struct
%
% ADDITIONALLY: we require the reference grid yRef,
%               use regularizer(...,'yRef',yRef)
%
% Output:
%    Sc    - hyperelastic energy
%    dS    - first derivative
%    d2S   - approximated Hessian
%
% see also hyperElastic.m
%==============================================================================

function [Sc,dS,d2S] = hyperElasticFEM(uc,yRef,Mesh,varargin)

persistent A  alphaLengthOld
if ~exist('A','var'),        A        = []; end;
if ~exist('omegaOld','var'), omegaOld = []; end;
if ~exist('mOld','var'),     mOld = []; end;
if nargin == 0, help(mfilename); runMinimalExample; return; end;

dS = []; d2S = [];

alpha       = regularizer('get','alpha');
if isempty(alpha), alpha=1; end
alphaLength = regularizer('get','alphaLength');
if isempty(alphaLength), alphaLength=1; end
alphaArea   = regularizer('get','alphaArea');
if isempty(alphaArea), alphaArea=1;     end
alphaVolume = regularizer('get','alphaVolume');
if isempty(alphaVolume), alphaVolume=1; end;
matrixFree  = 0;

doDerivative =( nargout>1);

for k=1:2:length(varargin), % overwrites default parameter
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

alphaLength= alpha*alphaLength;
alphaVolume= alpha*alphaVolume;

yc = yRef + uc;
dim = Mesh.dim;
vol = Mesh.vol;

if not(matrixFree), % matrix-based
    % =====================================================================
    % length regularization
    build = isempty(alphaLengthOld) || any(alphaLengthOld~=alphaLength) ||...
        isempty(A) || size(A,2)~=numel(uc);
    if build,
        alphaLengthOld = alphaLength;
        A = alphaLength * (Mesh.B'*sdiag(repmat(vol,[dim^2,1]))*Mesh.B);
    end
    dSlength  = uc'*A;
    Slength  = 0.5*dSlength*uc;
    d2Slength = A;
    
    GRAD = Mesh.GRAD;
    By = reshape(GRAD*reshape(yc,[],dim),[],dim^2);
    % =====================================================================
    % area regularization
    if dim==3,
        [cof, dCof] = cofactor3D(By,Mesh,doDerivative,matrixFree);
        % compute areas
        area = [
            sum(cof(:,[1 4 7]).^2,2);
            sum(cof(:,[2 5 8]).^2,2);
            sum(cof(:,[3 6 9]).^2,2);
            ];
        % compute penalty
        [H,dH,d2H] = phiDW(area,doDerivative);
        % compute area regularizer
        Sarea   = alphaArea * repmat(vol,[3 1])'*H;
        
        % derivative of area regularizer
        if doDerivative,
            dSarea = zeros(1,numel(yc));
            d2Sarea = sparse(numel(yc),numel(yc));
            dH = sdiag(vol)*reshape(dH,[],3); 
            d2H = reshape(d2H,[],3);
            for i=1:3,
               dA = 2*(sdiag(cof(:,i))*dCof{i}+sdiag(cof(:,i+3))*dCof{i+3}+sdiag(cof(:,i+6))*dCof{i+6});
               dSarea  = dSarea +   alphaArea*(dH(:,i)'*dA);
               d2Sarea = d2Sarea + dA'*sdiag(alphaArea*vol.*d2H(:,i))*dA;
            end
            clear dA dH d2H
%             dCof{4:9} = [];
        end
        
    else
        Sarea = 0;
        dSarea = 0;
        d2Sarea = 0;
    end
    
    % =====================================================================
    % volume regularization
    if dim==2,
        D1 = Mesh.dx1; D2 = Mesh.dx2;
        det = By(:,1).*By(:,4) - By(:,3) .* By(:,2);
        if doDerivative,
            dDet = [sdiag(By(:,4))*D1 + sdiag(-By(:,3))*D2, ...
                    sdiag(-By(:,2))*D1+ sdiag(By(:,1))*D2];
        end
    else
        D1 = Mesh.dx1; D2 = Mesh.dx2; D3 = Mesh.dx3;
        % built via cofactor (save some time)
        det = By(:,1).*cof(:,1)+By(:,2).*cof(:,2)+By(:,3).*cof(:,3);
       
        if doDerivative,
            Z = sparse(size(D1,1),size(D1,2));
            % simple product rule
            dDet = [sdiag(cof(:,1))*D1 + sdiag(cof(:,2))*D2 + sdiag(cof(:,3))*D3,Z,Z]...
                + sdiag(By(:,1))*dCof{1} + sdiag(By(:,2))*dCof{2} + sdiag(By(:,3))*dCof{3};
        end
    end
    [G,dG,d2G] = psi(det,doDerivative);
    Svolume   = alphaVolume *(vol'*G);
    Sc = Slength +Sarea + Svolume;

    if doDerivative
        % derivative of volume
        dSvolume  = alphaVolume * ((vol.*dG)'* dDet);
        d2Svolume = alphaVolume * (dDet' * (sdiag(vol.*d2G)) * dDet);
        if dim==3
            dS  = dSlength  + dSarea + dSvolume;
            d2S = d2Slength + d2Sarea + d2Svolume;
        elseif dim==2
            dS = dSlength + dSvolume;
            d2S = d2Slength + d2Svolume;
        end
    end
    
else % matrix-free
    d2S.regularizer = regularizer;
    d2S.alpha       = alpha;
    d2S.yc          = yc;
    d2S.solver      = 'PCG-hyperElastic';%@FEMMultiGridSolveHyper;
    
    % code only diffusion part for B
    d2S.By     = @(u,Mesh) Mesh.mfGRAD.D(u);
    d2S.BTy    = @(u,Mesh) Mesh.mfGRAD.Dadj(u);
    d2S.B      = @(Mesh)   Mesh.GRAD;
    
    % give seperate handles for diagonals of length, area and volume
    d2S.diagLength   = @(Mesh) getDiagLength(Mesh,alphaLength);
    d2S.diagArea     = [];%@(Mesh) getDiagArea(Mesh,alphaLength);
    d2S.diagVol      = @(Mesh) getDiagVolume(Mesh,yc,alphaVolume);
    d2S.diag         = @(yc) d2S.diagLength(Mesh) + d2S.diagVol(Mesh);
    
    % length regularizer
    dSlength = transpose( alphaLength* d2S.BTy(repmat(Mesh.vol,[dim^2,1]).*d2S.By(uc,Mesh),Mesh) );
    Slength  = .5*dSlength*uc;
    
    if doDerivative
        d2Slength = @(uc) alphaLength * d2S.BTy(repmat(Mesh.vol,[dim^2,1]).*d2S.By(uc,Mesh),Mesh);
    end
    
    % area regularizer
    if dim==3
        yc = reshape(yc,[],dim);
        [cof, dCof] = cofactor3D([],Mesh,doDerivative,matrixFree);
        % compute areas
        area = [
                cof{1}(yc).^2+cof{4}(yc).^2+cof{7}(yc).^2;
                cof{2}(yc).^2+cof{5}(yc).^2+cof{8}(yc).^2;
                cof{3}(yc).^2+cof{6}(yc).^2+cof{9}(yc).^2;
                ];    
        % compute penalty
        [H,dH,d2H] = phiDW(area,doDerivative);
        % compute area regularizer
        Sarea   = alphaArea * repmat(vol,[3 1])'*H;
        
        if doDerivative
            dH = reshape(dH,[],3);
            d2H = reshape(d2H,[],3);
            dAadj = @(x,i) 2*(dCof.dCofadj{i}(x.*cof{i}(yc),yc) + dCof.dCofadj{i+3}(x.*cof{i+3}(yc),yc) + dCof.dCofadj{i+6}(x.*cof{i+6}(yc),yc));
            dSarea = zeros(1,Mesh.nnodes*dim);
            for i=1:3          
               dSarea  = dSarea +   alphaArea*dAadj(vol.*dH(:,i),i)';
            end
            
            dA = @(x,i) 2*(cof{i}(yc).*dCof.dCof{i}(x,yc) + cof{i+3}(yc).*dCof.dCof{i+3}(x,yc) + cof{i+6}(yc).*dCof.dCof{i+6}(x,yc));
            d2Sarea = @(x) dAadj(alphaArea*vol.*d2H(:,1).*dA(reshape(x,[],3),1),1) + dAadj(alphaArea*vol.*d2H(:,2).*dA(reshape(x,[],3),2),2) + dAadj(alphaArea*vol.*d2H(:,3).*dA(reshape(x,[],3),3),3);
            
        end
        
    else
        Sarea = 0;
        dSarea = 0;
        d2Sarea = 0;
    end
    
    % volume regularizer
    if dim==2
        % volume term
        dx1 = Mesh.mfdx1; dx2 = Mesh.mfdx2;
        yc = reshape(yc,[],2);
        By = [dx1.D(yc(:,1))  dx2.D(yc(:,1)) dx1.D(yc(:,2)) dx2.D(yc(:,2))];
        det = By(:,1).*By(:,4) - By(:,3) .* By(:,2);
        [G,dG,d2G] = psi(det,doDerivative);
        Svolume = alphaVolume*sum(Mesh.vol.*G);
    
        if doDerivative
            
            dDet  = @(x)   By(:,4).*dx1.D(x(:,1))- By(:,3).*dx2.D(x(:,1))...
                                -By(:,2).*dx1.D(x(:,2))+By(:,1).*dx2.D(x(:,2));

            dDetadj = @(x) [ ...
                             dx1.Dadj(By(:,4).*x)-dx2.Dadj(By(:,3).*x);...
                            -dx1.Dadj(By(:,2).*x)+dx2.Dadj(By(:,1).*x)]; 
                        
            dSvolume = alphaVolume * transpose(dDetadj(dG.*Mesh.vol));
            
            d2Svolume = @(x) dDetadj(vol.*d2G.*dDet(reshape(x,[],2)));
            
        end
        
    else
        dx1 = Mesh.mfdx1; dx2 = Mesh.mfdx2; dx3 = Mesh.mfdx3;
        det = dx1.D(yc(:,1)).*cof{1}(yc) + dx2.D(yc(:,1)).*cof{2}(yc) + dx3.D(yc(:,1)).*cof{3}(yc);
        [G,dG,d2G] = psi(det,doDerivative);
        Svolume = alphaVolume*sum(Mesh.vol.*G);
        
        if doDerivative
            
            dDet = @(x)  cof{1}(yc).*dx1.D(x(:,1)) + cof{2}(yc).*dx2.D(x(:,1)) + cof{3}(yc).*dx3.D(x(:,1)) +...
                                  cof{4}(yc).*dx1.D(x(:,2)) + cof{5}(yc).*dx2.D(x(:,2)) + cof{6}(yc).*dx3.D(x(:,2)) +...
                                  cof{7}(yc).*dx1.D(x(:,3)) + cof{8}(yc).*dx2.D(x(:,3)) + cof{9}(yc).*dx3.D(x(:,3));
                              
            dDetadj = @(x) [ ...
                                 dx1.Dadj(cof{1}(yc).*x) + dx2.Dadj(cof{2}(yc).*x) + dx3.Dadj(cof{3}(yc).*x);...
                                 dx1.Dadj(cof{4}(yc).*x) + dx2.Dadj(cof{5}(yc).*x) + dx3.Dadj(cof{6}(yc).*x);...
                                 dx1.Dadj(cof{7}(yc).*x) + dx2.Dadj(cof{8}(yc).*x) + dx3.Dadj(cof{9}(yc).*x)
                                 ]; 
            
            dSvolume = alphaVolume * transpose(dDetadj(dG.*Mesh.vol));
            
            d2Svolume = @(x) alphaVolume * dDetadj(vol.*d2G.*dDet(reshape(x,[],3)));
            
        end
        
    end  
    Sc = Slength +Sarea + Svolume;
    
    if doDerivative
        if dim==3
            dS  = dSlength  + dSarea + dSvolume;
            d2S.d2S = @(x) d2Slength(x)  + d2Sarea(x) + d2Svolume(x);
        elseif dim==2
            dS = dSlength + dSvolume;
            d2S.d2S  = @(x) d2Slength(x) + d2Svolume(x);
        end
    end
    
    
end

function [cof,dCof] = cofactor3D(By,Mesh,doDerivative,matrixFree)
     
if ~matrixFree
    
    dCof = [];
    D = @(i,j) By(:,(j-1)*3+i); % gets partial_i y^j
    cof = [
        D(2,2).*D(3,3)-D(3,2).*D(2,3),...%ok
       -D(1,2).*D(3,3)+D(3,2).*D(1,3),...%trouble
        D(1,2).*D(2,3)-D(2,2).*D(1,3),...%ok
       -D(2,1).*D(3,3)+D(3,1).*D(2,3),...%ok
        D(1,1).*D(3,3)-D(3,1).*D(1,3),...%ok
       -D(1,1).*D(2,3)+D(2,1).*D(1,3),...%ok
        D(2,1).*D(3,2)-D(3,1).*D(2,2),...%ok
       -D(1,1).*D(3,2)+D(3,1).*D(1,2),...%ok
        D(1,1).*D(2,2)-D(2,1).*D(1,2)... %ok
        ];

    if doDerivative
        D1 = Mesh.dx1; D2 = Mesh.dx2; D3 = Mesh.dx3;
        Z = sparse(size(D1,1), size(D2,2));
        D = @(i,j) sdiag(By(:,(j-1)*3+i));

        dCof{1} = [  Z                 , D(3,3)*D2-D(2,3)*D3, D(2,2)*D3-D(3,2)*D2];
        dCof{2} = [  Z                 , D(1,3)*D3-D(3,3)*D1, D(3,2)*D1-D(1,2)*D3];
        dCof{3} = [  Z                 , D(2,3)*D1-D(1,3)*D2, D(1,2)*D2-D(2,2)*D1];
        dCof{4} = [ D(2,3)*D3-D(3,3)*D2, Z                  , D(3,1)*D2-D(2,1)*D3];
        dCof{5} = [ D(3,3)*D1-D(1,3)*D3, Z                  , D(1,1)*D3-D(3,1)*D1];
        dCof{6} = [ D(1,3)*D2-D(2,3)*D1, Z                  , D(2,1)*D1-D(1,1)*D2];
        dCof{7} = [ D(3,2)*D2-D(2,2)*D3 , D(2,1)*D3-D(3,1)*D2, Z];
        dCof{8} = [ D(1,2)*D3-D(3,2)*D1 , D(3,1)*D1-D(1,1)*D3, Z];
        dCof{9} = [ D(2,2)*D1-D(1,2)*D2 , D(1,1)*D2-D(2,1)*D1, Z];
    end

else
      
    dx1 = Mesh.mfdx1; dx2 = Mesh.mfdx2;  dx3 = Mesh.mfdx3;
            
    cof{1} = @(yc)  dx2.D(yc(:,2)).*dx3.D(yc(:,3)) - dx3.D(yc(:,2)).*dx2.D(yc(:,3));
    cof{2} = @(yc) -dx1.D(yc(:,2)).*dx3.D(yc(:,3)) + dx3.D(yc(:,2)).*dx1.D(yc(:,3));
    cof{3} = @(yc)  dx1.D(yc(:,2)).*dx2.D(yc(:,3)) - dx2.D(yc(:,2)).*dx1.D(yc(:,3));
    cof{4} = @(yc) -dx2.D(yc(:,1)).*dx3.D(yc(:,3)) + dx3.D(yc(:,1)).*dx2.D(yc(:,3));
    cof{5} = @(yc)  dx1.D(yc(:,1)).*dx3.D(yc(:,3)) - dx3.D(yc(:,1)).*dx1.D(yc(:,3));
    cof{6} = @(yc) -dx1.D(yc(:,1)).*dx2.D(yc(:,3)) + dx2.D(yc(:,1)).*dx1.D(yc(:,3));
    cof{7} = @(yc)  dx2.D(yc(:,1)).*dx3.D(yc(:,2)) - dx3.D(yc(:,1)).*dx2.D(yc(:,2));
    cof{8} = @(yc) -dx1.D(yc(:,1)).*dx3.D(yc(:,2)) + dx3.D(yc(:,1)).*dx1.D(yc(:,2));
    cof{9} = @(yc)  dx1.D(yc(:,1)).*dx2.D(yc(:,2)) - dx2.D(yc(:,1)).*dx1.D(yc(:,2));
    
    if doDerivative

        dx{1}.D = Mesh.mfdx1.D; dx{1}.Dadj = Mesh.mfdx1.Dadj;
        dx{2}.D = Mesh.mfdx2.D; dx{2}.Dadj = Mesh.mfdx2.Dadj;
        dx{3}.D = Mesh.mfdx3.D; dx{3}.Dadj = Mesh.mfdx3.Dadj;

        % adjoint/transpose of an element of dCof matrix
        dCofadjele = @(i,j,k,x,yc) dx{i}.Dadj(dx{j}.D(yc(:,k)).*x);

		dCof.dCofadj{1} = @(x,yc) [zeros(Mesh.nnodes,1); dCofadjele(2,3,3,x,yc)-dCofadjele(3,2,3,x,yc); dCofadjele(3,2,2,x,yc)-dCofadjele(2,3,2,x,yc)];
		dCof.dCofadj{2} = @(x,yc) [zeros(Mesh.nnodes,1); dCofadjele(3,1,3,x,yc)-dCofadjele(1,3,3,x,yc); dCofadjele(1,3,2,x,yc)-dCofadjele(3,1,2,x,yc)];
		dCof.dCofadj{3} = @(x,yc) [zeros(Mesh.nnodes,1); dCofadjele(1,2,3,x,yc)-dCofadjele(2,1,3,x,yc); dCofadjele(2,1,2,x,yc)-dCofadjele(1,2,2,x,yc)];

		dCof.dCofadj{4} = @(x,yc) [dCofadjele(3,2,3,x,yc)-dCofadjele(2,3,3,x,yc); zeros(Mesh.nnodes,1); dCofadjele(2,3,1,x,yc)-dCofadjele(3,2,1,x,yc)];
		dCof.dCofadj{5} = @(x,yc) [dCofadjele(1,3,3,x,yc)-dCofadjele(3,1,3,x,yc); zeros(Mesh.nnodes,1); dCofadjele(3,1,1,x,yc)-dCofadjele(1,3,1,x,yc)];
		dCof.dCofadj{6} = @(x,yc) [dCofadjele(2,1,3,x,yc)-dCofadjele(1,2,3,x,yc); zeros(Mesh.nnodes,1); dCofadjele(1,2,1,x,yc)-dCofadjele(2,1,1,x,yc)];

		dCof.dCofadj{7} = @(x,yc) [dCofadjele(2,3,2,x,yc)-dCofadjele(3,2,2,x,yc); dCofadjele(3,2,1,x,yc)-dCofadjele(2,3,1,x,yc); zeros(Mesh.nnodes,1)];
		dCof.dCofadj{8} = @(x,yc) [dCofadjele(3,1,2,x,yc)-dCofadjele(1,3,2,x,yc); dCofadjele(1,3,1,x,yc)-dCofadjele(3,1,1,x,yc); zeros(Mesh.nnodes,1)];
		dCof.dCofadj{9} = @(x,yc) [dCofadjele(1,2,2,x,yc)-dCofadjele(2,1,2,x,yc); dCofadjele(2,1,1,x,yc)-dCofadjele(1,2,1,x,yc); zeros(Mesh.nnodes,1)];

		% element of dCof matrix
		dCofele = @(i,j,k,l,x,yc) dx{i}.D(yc(:,j)).*dx{k}.D(x(:,l));

		dCof.dCof{1} = @(x,yc) dCofele(3,3,2,2,x,yc)-dCofele(2,3,3,2,x,yc) + dCofele(2,2,3,3,x,yc) - dCofele(3,2,2,3,x,yc);
		dCof.dCof{2} = @(x,yc) dCofele(1,3,3,2,x,yc)-dCofele(3,3,1,2,x,yc) + dCofele(3,2,1,3,x,yc) - dCofele(1,2,3,3,x,yc);
		dCof.dCof{3} = @(x,yc) dCofele(2,3,1,2,x,yc)-dCofele(1,3,2,2,x,yc) + dCofele(1,2,2,3,x,yc) - dCofele(2,2,1,3,x,yc);

		dCof.dCof{4} = @(x,yc) dCofele(2,3,3,1,x,yc)-dCofele(3,3,2,1,x,yc) + dCofele(3,1,2,3,x,yc) - dCofele(2,1,3,3,x,yc);
		dCof.dCof{5} = @(x,yc) dCofele(3,3,1,1,x,yc)-dCofele(1,3,3,1,x,yc) + dCofele(1,1,3,3,x,yc) - dCofele(3,1,1,3,x,yc);
		dCof.dCof{6} = @(x,yc) dCofele(1,3,2,1,x,yc)-dCofele(2,3,1,1,x,yc) + dCofele(2,1,1,3,x,yc) - dCofele(1,1,2,3,x,yc);

		dCof.dCof{7} = @(x,yc) dCofele(3,2,2,1,x,yc)-dCofele(2,2,3,1,x,yc) + dCofele(2,1,3,2,x,yc) - dCofele(3,1,2,2,x,yc);
		dCof.dCof{8} = @(x,yc) dCofele(1,2,3,1,x,yc)-dCofele(3,2,1,1,x,yc) + dCofele(3,1,1,2,x,yc) - dCofele(1,1,3,2,x,yc);
		dCof.dCof{9} = @(x,yc) dCofele(2,2,1,1,x,yc)-dCofele(1,2,2,1,x,yc) + dCofele(1,1,2,2,x,yc) - dCofele(2,1,1,2,x,yc);
       
    end
    
end

function D = getDiagVolume(Mesh,yc,alphaVolume)
dim = Mesh.dim;
vol = Mesh.vol;
if dim==2
    xn = Mesh.xn;
    % get edges
    e1 =  Mesh.mfPi(xn,1) - Mesh.mfPi(xn,3);
    e2 =  Mesh.mfPi(xn,2) - Mesh.mfPi(xn,3);
    e3 =  Mesh.mfPi(xn,1) - Mesh.mfPi(xn,2);
    
    % compute gradients of basis functions
    dphi1 =   ([ e2(:,2) -e2(:,1)]./[2*vol,2*vol]);
    dphi2 =   ([-e1(:,2)  e1(:,1)]./[2*vol,2*vol]);
    dphi3 =   ([ e3(:,2) -e3(:,1)]./[2*vol,2*vol]);
    
    % compute diagonal of volume regularizer
    % MB : dDet = [sdiag(By(:,4)), sdiag(-By(:,3)) , sdiag(-By(:,2)), sdiag(By(:,1))] * B;
    [dx1, dx2] = getGradientMatrixFEM(Mesh,1);
    yc = reshape(yc,[],2);
    By = [dx1.D(yc(:,1))  dx2.D(yc(:,1)) dx1.D(yc(:,2)) dx2.D(yc(:,2))];
    det = By(:,1).*By(:,4) - By(:,3) .* By(:,2);
    [~,~,d2G] = psi(det,1);
    
    % get boundaries
    Dxi = @(i,j) Mesh.mfPi(vol.*d2G.*(By(:,j).*dphi1(:,i)).^2,1) ... 
        + Mesh.mfPi(vol.*d2G.*(By(:,j).*dphi2(:,i)).^2,2) ...
        + Mesh.mfPi(vol.*d2G.*(By(:,j).*dphi3(:,i)).^2,3); % diagonal of Dxi'*Dxi 
    
    Dxy = @(i,j,k,l) Mesh.mfPi(vol.*d2G.*(By(:,j).*dphi1(:,i)).*(By(:,l).*dphi1(:,k)),1) ... % byproduct terms for verifications
        + Mesh.mfPi(vol.*d2G.*(By(:,j).*dphi2(:,i)).*(By(:,l).*dphi2(:,k)),2) ...
        + Mesh.mfPi(vol.*d2G.*(By(:,j).*dphi3(:,i)).*(By(:,l).*dphi3(:,k)),3); % diagonal of Dxi'*Dxi 
    
    D = [Dxi(1,4)+ Dxi(2,3) - 2*Dxy(1,4,2,3) ;Dxi(1,2)+ Dxi(2,1) - 2*Dxy(1,2,2,1)];
else
    error('nyi');
end
D = alphaVolume*D;


% compute d2Svol*x
function By = volumeOperator(yc,Mesh,x,vol)
[dx1, dx2] = getGradientMatrixFEM(Mesh,1);
% volume regularization
yc = reshape(yc,[],2);
By = [dx1.D(yc(:,1))  dx2.D(yc(:,1)) dx1.D(yc(:,2)) dx2.D(yc(:,2))];
det = By(:,1).*By(:,4) - By(:,3) .* By(:,2);
[~,~,d2G] = psi(det,1);

%      MB : dDet = [sdiag(By(:,4)), sdiag(-By(:,3)) , sdiag(-By(:,2)), sdiag(By(:,1))] * B;

dDet  = @(x)   By(:,4).*dx1.D(x(:,1))- By(:,3).*dx2.D(x(:,1))...
    -By(:,2).*dx1.D(x(:,2))+By(:,1).*dx2.D(x(:,2));

dDetAdj  = @(x) [ ...
    dx1.Dadj(By(:,4).*x)-dx2.Dadj(By(:,3).*x);...
    -dx1.Dadj(By(:,2).*x)+dx2.Dadj(By(:,1).*x)
    ];
By = dDetAdj(vol.*d2G.*dDet(reshape(x,[],2)));



function D = getDiagLength(Mesh,alphaLength)
dim = Mesh.dim;
vol = Mesh.vol;
if dim==2
    xn = Mesh.xn;
    % get edges
    e1 =  Mesh.mfPi(xn,1) - Mesh.mfPi(xn,3);
    e2 =  Mesh.mfPi(xn,2) - Mesh.mfPi(xn,3);
    e3 =  Mesh.mfPi(xn,1) - Mesh.mfPi(xn,2);
    
    % compute gradients of basis functions
    dphi1 =   ([ e2(:,2) -e2(:,1)]./[2*vol,2*vol]);
    dphi2 =   ([-e1(:,2)  e1(:,1)]./[2*vol,2*vol]);
    dphi3 =   ([ e3(:,2) -e3(:,1)]./[2*vol,2*vol]);
    
    % get boundaries
    Dxi = @(i) Mesh.mfPi(vol.*dphi1(:,i).^2,1) + Mesh.mfPi(vol.*dphi2(:,i).^2,2) + Mesh.mfPi(vol.*dphi3(:,i).^2,3); % diagonal of Dx1'*Dx1
        
    D = Dxi(1)+ Dxi(2);
    D = [D;D];
else
    error('nyi');
end
D = alphaLength*D(:);



function [G dG d2G] = psi(x,doDerivative)
%
% psi(x) = ((x-1)^2/x)^2
%
% psi satisfies the three important conditions
%      psi(x) > 0, forall x
%      psi(1) = 0
%      psi is convex
%      psi(x) = psi(1/x)
%      psi yields det(Dy) in L_2
dG = [];
d2G = [];

G = (x-1).*(x-1) ./x;
G = G.*G;
if doDerivative,
    dG  = 2* (x-1).^3 .* (x+1)./ x.^3;
    d2G = 2* (x.^4-4*x+3) ./ x.^4;
end

function [G dG d2G] = phiC(x,doDerivative)
%
%   phiC is a convex penalty function for surface area
%   regularization. This function is needed in the existence proof. In
%   order to be convex, only area growth can be penalized
%
%
%    phi(x>=1) = 0.5 * (A/ ARef - 1)^2
%    phi(x<1 ) = 0
%
%    A    - area after deformation (24 Triangles per voxel, scalar)
%    ARef - area of reference configuration (24 Triangles per voxel, scalar)
%
dG = [];
d2G = [];
G = 0.5* (x-1).*(x-1);
G(x<1)=0;
if doDerivative,
    dG = (x-1);
    dG(x<1)=0;
    d2G = ones(size(dG));
    d2G(x<1)=0;
end

function [G dG d2G] = phiDW(x,doDerivative)
%
%   phiDW is a penalty function for surface area
%   regularization that penalizes growth and shrinkage of area. However,
%   this function is a double well and thus not convex.
%
%    phi(A) = 0.5 * (A/ ARef - 1)^2
%
%    A    - area after deformation (24 Triangles per voxel, scalar)
%    ARef - area of reference configuration (24 Triangles per voxel, scalar)
dG = [];
d2G = [];
G = 0.5* (x-1).*(x-1);
if doDerivative,
    dG = (x-1);
    d2G = ones(size(dG));
end

% shortcut for spdiags
function D = sdiag(v)
D = spdiags(v(:),0,numel(v),numel(v));

function runMinimalExample
%========= 2D
omega = [0,10,0,8]; m = [17,16]; p = [5,6];
w = zeros([p,2]);  w(3,3,1) = 0.06; w(3,4,2) = -0.05;
Mesh   = TriMesh1(omega,m);
xn = Mesh.xn;
yn = splineTransformation2D(w(:),xn(:),'omega',omega,'m',m+1,'p',p,'Q',[]);

regOptn = {'alpha',1,'alphaLength',0,'alphaArea',0,'alphaVolume',1};

fctn = @(y) feval(mfilename,y-xn(:),xn(:),Mesh,regOptn{:},'matrixFree',0);
[Sc, dS, d2S] = fctn(yn(:));
checkDerivative(fctn,yn(:));

%========== 3D
omega = [0,10,0,8,0 4]; m = [5,6,8]; type = 1;
Mesh = TetraMesh1(omega,m);
xn = Mesh.xn;
yn = xn + 1e-1* randn(size(xn));

regOptn = {'alpha',1,'alphaLength',0,'alphaArea',1,'alphaVolume',0};

fctn = @(y) feval(mfilename,y-xn(:),xn(:),Mesh,regOptn{:},'matrixFree',0);
[Sc, dS, d2S] = fctn(yn(:));
checkDerivative(fctn,yn(:));

fctn = @(y) feval(mfilename,y-xn(:),xn(:),Mesh,regOptn{:},'matrixFree',1);
[ScMF, dSMF, d2SMF] = fctn(yn(:));


