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
            %Z = sparse(size(D1,1),size(D1,2));
            % simple product rule
            %dDet = [sdiag(cof(:,1))*D1 + sdiag(cof(:,2))*D2 + sdiag(cof(:,3))*D3,Z,Z]...
            %    + sdiag(By(:,1))*dCof{1} + sdiag(By(:,2))*dCof{2} + sdiag(By(:,3))*dCof{3};

            dDet = [sdiag(cof(:,1))*D1 + sdiag(cof(:,2))*D2 + sdiag(cof(:,3))*D3,...
                    sdiag(cof(:,4))*D1 + sdiag(cof(:,5))*D2 + sdiag(cof(:,6))*D3,...
                    sdiag(cof(:,7))*D1 + sdiag(cof(:,8))*D2 + sdiag(cof(:,9))*D3];

        end
    end
    [G,dG,d2G] = psi(det,doDerivative);
    Svolume   = alphaVolume *(vol'*G);
    Sc = Slength + Sarea + Svolume;

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
    By         = reshape(Mesh.mfGRAD.D(yc),[],dim^2); % Compute once and re-use (function calls too costly)
    
    
    %======================================================================
    % Length regularizer
    dSlength = transpose( alphaLength* d2S.BTy(repmat(Mesh.vol,[dim^2,1]).*d2S.By(uc,Mesh),Mesh) );
    Slength  = .5*dSlength*uc;
    
    if doDerivative
        d2Slength = @(uc) alphaLength * d2S.BTy(repmat(Mesh.vol,[dim^2,1]).*d2S.By(uc,Mesh),Mesh);
        d2SlengthDiag   = @(Mesh) getDiagLength(Mesh,alphaLength);
    end
    
    %======================================================================
    % Area regularizer
    if dim==3
        yc = reshape(yc,[],dim);
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
        
        if doDerivative
            dH = reshape(dH,[],3);
            d2H = reshape(d2H,[],3);
            dAadj   = @(x,i) 2*(dCof.dCofadj{i}(x.*cof(:,i),yc) + dCof.dCofadj{i+3}(x.*cof(:,i+3),yc) + dCof.dCofadj{i+6}(x.*cof(:,i+6),yc));
            dSarea = zeros(1,Mesh.nnodes*dim);
            for i=1:3          
               dSarea  = dSarea +   alphaArea*dAadj(vol.*dH(:,i),i)';
            end
            
            dA          = @(x,i) 2*(cof(:,i).*dCof.dCof{i}(x,yc) + cof(:,i+3).*dCof.dCof{i+3}(x,yc) + cof(:,i+6).*dCof.dCof{i+6}(x,yc));
            d2Sarea     = @(x) dAadj(alphaArea*vol.*d2H(:,1).*dA(reshape(x,[],3),1),1) + dAadj(alphaArea*vol.*d2H(:,2).*dA(reshape(x,[],3),2),2)...
                             + dAadj(alphaArea*vol.*d2H(:,3).*dA(reshape(x,[],3),3),3);  
            d2SareaDiag = @(Mesh) getDiagArea(Mesh, alphaArea, By, cof);
        end        
    else
        Sarea = 0;
        dSarea = 0;
        d2Sarea = 0;
    end
    
    %======================================================================
    % Volume regularizer
    if dim==2
        % volume term
        dx1 = Mesh.mfdx1; dx2 = Mesh.mfdx2;
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
            d2SvolumeDiag   = @(Mesh) getDiagVolume(Mesh,alphaVolume, By, det, []);
        end   
    else % (dim==3)
        dx1 = Mesh.mfdx1; dx2 = Mesh.mfdx2; dx3 = Mesh.mfdx3;
        det     = By(:,1).*cof(:,1) + By(:,2).*cof(:,2) + By(:,3).*cof(:,3);
        [G,dG,d2G] = psi(det,doDerivative);
        Svolume = alphaVolume*sum(Mesh.vol.*G);
        
        if doDerivative
            dDet    = @(x)  cof(:,1).*dx1.D(x(:,1)) + cof(:,2).*dx2.D(x(:,1)) + cof(:,3).*dx3.D(x(:,1)) ...
                          + cof(:,4).*dx1.D(x(:,2)) + cof(:,5).*dx2.D(x(:,2)) + cof(:,6).*dx3.D(x(:,2)) ...
                          + cof(:,7).*dx1.D(x(:,3)) + cof(:,8).*dx2.D(x(:,3)) + cof(:,9).*dx3.D(x(:,3));
                              
            dDetadj = @(x) [ ...
                               dx1.Dadj(cof(:,1).*x) + dx2.Dadj(cof(:,2).*x) + dx3.Dadj(cof(:,3).*x);...
                               dx1.Dadj(cof(:,4).*x) + dx2.Dadj(cof(:,5).*x) + dx3.Dadj(cof(:,6).*x);...
                               dx1.Dadj(cof(:,7).*x) + dx2.Dadj(cof(:,8).*x) + dx3.Dadj(cof(:,9).*x)...
                           ]; 
            
            dSvolume = alphaVolume * transpose(dDetadj(dG.*Mesh.vol));
            
            d2Svolume = @(x) alphaVolume * dDetadj(vol.*d2G.*dDet(reshape(x,[],3)));
            d2SvolumeDiag   = @(Mesh) getDiagVolume(Mesh,alphaVolume, By, det, cof);
        end
        
    end  
    Sc = Slength +Sarea + Svolume;
    
    if doDerivative
        if dim==3
            dS  = dSlength  + dSarea + dSvolume;
            d2S.d2S = @(x) d2Slength(x)  + d2Sarea(x) + d2Svolume(x);
            d2S.d2Slength       = @(x) d2Slength(x);
            d2S.d2Sarea         = @(x) d2Sarea(x);
            d2S.d2Svolume       = @(x) d2Svolume(x);
        
            % Give seperate handles for diagonals of length, area and volume
            d2S.d2Sdiag         = @(yc) d2SlengthDiag(Mesh) + d2SareaDiag(Mesh) + d2SvolumeDiag(Mesh);
            d2S.d2SlengthDiag   = @(yc) d2SlengthDiag(Mesh);
            d2S.d2SareaDiag     = @(yc) d2SareaDiag(Mesh);
            d2S.d2SvolumeDiag   = @(yc) d2SvolumeDiag(Mesh);

        elseif dim==2
            dS = dSlength + dSvolume;
            d2S.d2S  = @(x) d2Slength(x) + d2Svolume(x);
            d2S.d2Slength       = @(x) d2Slength(x);
            d2S.d2Svolume       = @(x) d2Svolume(x);
            
            % Give seperate handles for diagonals of length and volume
            d2S.d2Sdiag         = @(yc) d2SlengthDiag(Mesh) + d2SvolumeDiag(Mesh);
            d2S.d2SlengthDiag   = @(yc) d2SlengthDiag(Mesh);
            d2S.d2SvolumeDiag   = @(yc) d2SvolumeDiag(Mesh);
            
        end
    end
    
    
end
end

function D = getDiagLength(Mesh,alphaLength)
%
%   Computes diagonal of length Hessian for matrix-free implementation.  
%
%   Input:         Mesh - mesh structure
%           alphaLength - regularization parameter for length term
%
%   Output:         D - diagonal of length term
%

dim = Mesh.dim;
vol = Mesh.vol;
if dim==2
    
    dphi  = Mesh.dphi;    
    dphi1 = dphi{1}; dphi2 = dphi{2}; dphi3 = dphi{3};
    
    % get boundaries
    Dxi = @(i) Mesh.mfPi(vol.*dphi1(:,i).^2,1) + Mesh.mfPi(vol.*dphi2(:,i).^2,2) + Mesh.mfPi(vol.*dphi3(:,i).^2,3); % diagonal of Dx1'*Dx1
        
    D = Dxi(1)+ Dxi(2);
    D = [D;D];
else
    dphi  = Mesh.dphi;    
    dphi1 = dphi{1}; dphi2 = dphi{2}; dphi3 = dphi{3}; dphi4 = dphi{4};
    
    % get boundaries
    Dxi = @(i) Mesh.mfPi(vol.*dphi1(:,i).^2,1) + ...
                Mesh.mfPi(vol.*dphi2(:,i).^2,2) + ...
                    Mesh.mfPi(vol.*dphi3(:,i).^2,3) + ...
                        Mesh.mfPi(vol.*dphi4(:,i).^2,4); % diagonal of Dx1'*Dx1
        
    D = Dxi(1)+ Dxi(2) + Dxi(3);
    D = [D;D;D];

end
D = alphaLength*D(:);
end


function D = getDiagArea(Mesh,alphaArea,By,cof)
%
%   Computes diagonal of surface area Hessian for matrix-free
%   implementation. Only needed for dim==3 formulation 
%
%   Input:       Mesh - mesh structure
%           alphaArea - regularization parameter for area term
%                  By - #nodes x 9 matrix with partial derivatives of
%                           dimensions of y w.r.t. each direction
%                 cof - #tri x 9 matrix with pre-computed cofactors of y
%
%   Output:         D - diagonal of area term
%
dim = Mesh.dim;
vol = Mesh.vol;

dphi = Mesh.dphi;    
dphi1 = dphi{1}; dphi2 = dphi{2}; dphi3 = dphi{3}; dphi4 = dphi{4};

% Note: for double-well function can d2H = identity so could omit these
% terms. Necessary for phiC, so included here. 

% Compute areas
area = [
        sum(cof(:,[1 4 7]).^2,2);
        sum(cof(:,[2 5 8]).^2,2);
        sum(cof(:,[3 6 9]).^2,2);
       ];   
% Compute penalty
[~,~,d2H] = phiDW(area,1);
d2H = reshape(d2H,[],3); 

% Get boundaries.             
Dxi = @(i,j,k,l,m,y) Mesh.mfPi(vol.*d2H(:,y).*(cof(:,m).^2).*(By(:,j).*dphi1(:,i) - By(:,l).*dphi1(:,k)).^2,1) ... 
                 + Mesh.mfPi(vol.*d2H(:,y).*(cof(:,m).^2).*(By(:,j).*dphi2(:,i) - By(:,l).*dphi2(:,k)).^2,2) ... 
                 + Mesh.mfPi(vol.*d2H(:,y).*(cof(:,m).^2).*(By(:,j).*dphi3(:,i) - By(:,l).*dphi3(:,k)).^2,3) ... 
                 + Mesh.mfPi(vol.*d2H(:,y).*(cof(:,m).^2).*(By(:,j).*dphi4(:,i) - By(:,l).*dphi4(:,k)).^2,4); 

Dxy = @(i,j,k,l,m,y,p,q,r,s,n,z) Mesh.mfPi(vol.*d2H(:,y).*cof(:,m).*d2H(:,z).*cof(:,n).*(By(:,j).*dphi1(:,i) - By(:,l).*dphi1(:,k)).*(By(:,q).*dphi1(:,p) - By(:,s).*dphi1(:,r)),1) ...                
                           + Mesh.mfPi(vol.*d2H(:,y).*cof(:,m).*d2H(:,z).*cof(:,n).*(By(:,j).*dphi2(:,i) - By(:,l).*dphi2(:,k)).*(By(:,q).*dphi2(:,p) - By(:,s).*dphi2(:,r)),2) ...                
                           + Mesh.mfPi(vol.*d2H(:,y).*cof(:,m).*d2H(:,z).*cof(:,n).*(By(:,j).*dphi3(:,i) - By(:,l).*dphi3(:,k)).*(By(:,q).*dphi3(:,p) - By(:,s).*dphi3(:,r)),3) ...                
                           + Mesh.mfPi(vol.*d2H(:,y).*cof(:,m).*d2H(:,z).*cof(:,n).*(By(:,j).*dphi4(:,i) - By(:,l).*dphi4(:,k)).*(By(:,q).*dphi4(:,p) - By(:,s).*dphi4(:,r)),4);                
                

D1 = Dxi(3,8,2,9,4,1) + Dxi(1,9,3,7,5,2) + Dxi(2,7,1,8,6,3) + Dxi(2,6,3,5,7,1) + Dxi(3,4,1,6,8,2) + Dxi(1,5,2,4,9,3) ... %Square terms
   + 2*Dxy(3,8,2,9,4,1,2,6,3,5,7,1) + 2*Dxy(1,9,3,7,5,2,3,4,1,6,8,2) + 2*Dxy(2,7,1,8,6,3,1,5,2,4,9,3); %Cross terms

D2 = Dxi(2,9,3,8,1,1) + Dxi(3,7,1,9,2,2) + Dxi(1,8,2,7,3,3) + Dxi(3,2,2,3,7,1) + Dxi(1,3,3,1,8,2) + Dxi(2,1,1,2,9,3) ... %Square terms
   + 2*Dxy(2,9,3,8,1,1,3,2,2,3,7,1) + 2*Dxy(3,7,1,9,2,2,1,3,3,1,8,2) + 2*Dxy(1,8,2,7,3,3,2,1,1,2,9,3); %Cross terms

D3 = Dxi(3,5,2,6,1,1) + Dxi(1,6,3,4,2,2) + Dxi(2,4,1,5,3,3) + Dxi(2,3,3,2,4,1) + Dxi(3,1,1,3,5,2) + Dxi(1,2,2,1,6,3) ... %Square terms
   + 2*Dxy(3,5,2,6,1,1,2,3,3,2,4,1) + 2*Dxy(1,6,3,4,2,2,3,1,1,3,5,2) + 2*Dxy(2,4,1,5,3,3,1,2,2,1,6,3); %Cross terms
 
D = 4*alphaArea*[D1; D2; D3];

end


function D = getDiagVolume(Mesh,alphaVolume,By,det,cof)
%
%   Computes diagonal of volume Hessian for matrix-free
%   implementation.
%
%   Input:         Mesh - mesh structure
%           alphaVolume - regularization parameter for volume term
%                    By - #nodes x 9 matrix with partial derivatives of
%                           dimensions of y w.r.t. each direction
%                   det - #tri vector with determinant of mesh triangles/tetrahedra
%                   cof - #tri x 9 matrix with pre-computed cofactors of y
%
%   Output:         D - diagonal of volume term
%


dim = Mesh.dim;
vol = Mesh.vol;
if dim==2
    
    dphi  = Mesh.dphi;    
    dphi1 = dphi{1}; dphi2 = dphi{2}; dphi3 = dphi{3};
    
    % Compute diagonal term Psi"(detBy)
    [~,~,d2G] = psi(det,1);

    % get boundaries
    % Square terms
    Dxi = @(i,j) Mesh.mfPi(vol.*d2G.*(By(:,j).*dphi1(:,i)).^2,1) ... 
        + Mesh.mfPi(vol.*d2G.*(By(:,j).*dphi2(:,i)).^2,2) ...
        + Mesh.mfPi(vol.*d2G.*(By(:,j).*dphi3(:,i)).^2,3); % diagonal of Dxi'*Dxi 
    
    % Cross terms
    Dxy = @(i,j,k,l) Mesh.mfPi(vol.*d2G.*(By(:,j).*dphi1(:,i)).*(By(:,l).*dphi1(:,k)),1) ... % byproduct terms for verifications
        + Mesh.mfPi(vol.*d2G.*(By(:,j).*dphi2(:,i)).*(By(:,l).*dphi2(:,k)),2) ...
        + Mesh.mfPi(vol.*d2G.*(By(:,j).*dphi3(:,i)).*(By(:,l).*dphi3(:,k)),3); % diagonal of Dxi'*Dxi 
    
    D = [Dxi(1,4)+ Dxi(2,3) - 2*Dxy(1,4,2,3) ;Dxi(1,2)+ Dxi(2,1) - 2*Dxy(1,2,2,1)];

else
    dphi  = Mesh.dphi;   
    dphi1 = dphi{1}; dphi2 = dphi{2}; dphi3 = dphi{3}; dphi4 = dphi{4};
    
    % Compute diagonal term vol.*Psi"(detBy)
    [~,~,d2G] = psi(det,1);
    vd2G = vol.*d2G;
    
    % get boundaries
    % Square terms
    Dxi = @(i,j) Mesh.mfPi(vd2G.*(cof(:,j).*dphi1(:,i)).^2,1) ... 
               + Mesh.mfPi(vd2G.*(cof(:,j).*dphi2(:,i)).^2,2) ...
               + Mesh.mfPi(vd2G.*(cof(:,j).*dphi3(:,i)).^2,3) ... 
               + Mesh.mfPi(vd2G.*(cof(:,j).*dphi4(:,i)).^2,4); 
    % Cross terms
    Dxy = @(i,j,k,l)   Mesh.mfPi(vd2G.*(cof(:,j).*dphi1(:,i)).*(cof(:,l).*dphi1(:,k)),1) ... % byproduct terms for verifications
                     + Mesh.mfPi(vd2G.*(cof(:,j).*dphi2(:,i)).*(cof(:,l).*dphi2(:,k)),2) ...
                     + Mesh.mfPi(vd2G.*(cof(:,j).*dphi3(:,i)).*(cof(:,l).*dphi3(:,k)),3) ...
                     + Mesh.mfPi(vd2G.*(cof(:,j).*dphi4(:,i)).*(cof(:,l).*dphi4(:,k)),4); % diagonal of Dxi'*Dxi
       
    D = [Dxi(1,1) + Dxi(2,2) + Dxi(3,3) + 2*Dxy(1,1,2,2) + 2*Dxy(1,1,3,3) + 2*Dxy(2,2,3,3); ...
         Dxi(1,4) + Dxi(2,5) + Dxi(3,6) + 2*Dxy(1,4,2,5) + 2*Dxy(1,4,3,6) + 2*Dxy(2,5,3,6); ...
         Dxi(1,7) + Dxi(2,8) + Dxi(3,9) + 2*Dxy(1,7,2,8) + 2*Dxy(1,7,3,9) + 2*Dxy(2,8,3,9)];    
end
D = alphaVolume*D;
end


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
end

% shortcut for spdiags
function D = sdiag(v)
D = spdiags(v(:),0,numel(v),numel(v));
end

function runMinimalExample
%========= 2D
omega = [0,10,0,8]; m = [17,16]; p = [5,6];
w = zeros([p,2]);  w(3,3,1) = 0.06; w(3,4,2) = -0.05;
Mesh   = TriMesh1(omega,m);
xn = Mesh.xn;
yn = splineTransformation2D(w(:),xn(:),'omega',omega,'m',m+1,'p',p,'Q',[]);

regOptn = {'alpha',1,'alphaLength',1,'alphaArea',1,'alphaVolume',1};

fctn = @(y) feval(mfilename,y-xn(:),xn(:),Mesh,regOptn{:},'matrixFree',0);
[Sc, dS, d2S] = fctn(yn(:));
checkDerivative(fctn,yn(:));

fctn = @(y) feval(mfilename,y-xn(:),xn(:),Mesh,regOptn{:},'matrixFree',1);
[ScMF, dSMF, d2SMF] = fctn(yn(:));

% Check that matrix-based and matrix-free versions agree
p       = randn(size(yn(:)));
d2Sp    = d2S*p;
d2SMFp  = d2SMF.d2S(p);
d2Sdiagp    = diag(d2S);
d2SMFdiagp  = d2SMF.d2Sdiag(ones(size(p)));

fprintf('norm(d2Sp - d2SMFp) = %1.4e \n', norm(d2Sp(:) - d2SMFp(:)));
fprintf('norm(d2Sdiagp - d2SMFdiagp) = %1.4e \n', norm(d2Sdiagp(:) - d2SMFdiagp(:)));

%========== 3D
omega = [0,10,0,8,0 4]; m = [5,6,8]; type = 1;
%omega = [0,5,0,5,0 5]; m = [5,5,5]; type = 1;
Mesh = TetraMesh1(omega,m);
xn = Mesh.xn;
yn = xn + 1e-1* randn(size(xn));

regOptn = {'alpha',1,'alphaLength',1,'alphaArea',1,'alphaVolume',1};

fctn = @(y) feval(mfilename,y-xn(:),xn(:),Mesh,regOptn{:},'matrixFree',0);
[Sc, dS, d2S] = fctn(yn(:));
checkDerivative(fctn,yn(:));

fctn = @(y) feval(mfilename,y-xn(:),xn(:),Mesh,regOptn{:},'matrixFree',1);
[ScMF, dSMF, d2SMF] = fctn(yn(:));

% Check that matrix-based and matrix-free versions agree
p       = randn(size(yn(:)));
d2Sp    = d2S*p;
d2SMFp  = d2SMF.d2S(p);
d2Sdiagp    = diag(d2S);
d2SMFdiagp  = d2SMF.d2Sdiag(ones(size(p)));

fprintf('norm(d2Sp - d2SMFp) = %1.4e \n', norm(d2Sp(:) - d2SMFp(:)));
fprintf('norm(d2Sdiagp - d2SMFdiagp) = %1.4e \n', norm(d2Sdiagp(:) - d2SMFdiagp(:)));

% Check Hessian for 3D problems
h = logspace(-1,-10,10);
v = randn(size(yn(:))); 
f0   = Sc;
dvf  = reshape(dS,1,[])*v;
d2vf = v'*d2SMF.d2S(v);

for j=1:length(h)
  ft = feval(fctn,yn(:)+h(j)*v);      % function value
  T0(j) = norm(f0-ft);             % TaylorPoly 0
  T1(j) = norm(f0+h(j)*dvf - ft);  % TaylorPoly 1
  T2(j) = norm(f0+h(j)*dvf +0.5*h(j)*h(j)*d2vf - ft);  % TaylorPoly 1
  fprintf('h=%12.4e     T0=%12.4e    T1=%12.4e  T2=%12.4e \n',h(j),T0(j),T1(j),T2(j));
end

figure(); loglog(h,[T0;T1;T2]);
end



