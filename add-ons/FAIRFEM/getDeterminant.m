%==============================================================================
% This code is part of the Finite Element Method app for the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR/FAIRFEM 
%==============================================================================
%
% function [det,dDet] = getDeterminant(Mesh,yc,varargin)
%
% computes determinant of Jacobian of transformation discreteized on
% triangular or tetrahedral mesh.
%
% Input:
%   Mesh - description of mesh (struct)
%   yc   - nodes of finite element grid
%
% Output:
%   det  - Jacobian determinant for each element
%  dDet  - derivative w.r.t. yc
%
%==============================================================================
function [det,dDet] = getDeterminant(Mesh,yc,varargin)
dDet = [];
if nargin==0, help(mfilename); runMinimalExample; return; end;

matrixFree = regularizer('get','matrixFree');
dim = Mesh.dim;
yc = reshape(yc,[],dim);

doDerivative = (nargout>1);
for k=1:2:length(varargin)     % overwrites default parameter
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

switch dim
    case 2
        if ~matrixFree
            dx1 = Mesh.dx1; dx2 = Mesh.dx2;

            D1y = dx1*yc;
            D2y = dx2*yc;
            % compute volume
            det = D1y(:,1).*D2y(:,2) - D2y(:,1).*D1y(:,2);
        
            if doDerivative
                w = [D2y(:,2), D1y(:,2), D1y(:,1), D2y(:,1)];
                dDet =[
                    sdiag(w(:,1))*(dx1)-sdiag(w(:,2))*(dx2),...
                    sdiag(w(:,3))*(dx2)-sdiag(w(:,4))*(dx1)];
            end
        else
            dx1 = Mesh.mfdx1; dx2 = Mesh.mfdx2;
            By = [dx1.D(yc(:,1))  dx2.D(yc(:,1)) dx1.D(yc(:,2)) dx2.D(yc(:,2))];
            det = By(:,1).*By(:,4) - By(:,3) .* By(:,2);
            
            dDet.dDet = @(x) By(:,4).*dx1.D(x(:,1)) - By(:,3).*dx2.D(x(:,1)) + By(:,1).*dx2.D(x(:,2)) - By(:,2).*dx1.D(x(:,2));
            
            dDet.dDetadj  = @(x) [ ...
                dx1.Dadj(By(:,4).*x)-dx2.Dadj(By(:,3).*x);...
                -dx1.Dadj(By(:,2).*x)+dx2.Dadj(By(:,1).*x)]; 
        end
    case 3
        
        if ~matrixFree
            % compute derivative
            GRAD = Mesh.GRAD;
            By = reshape(GRAD*reshape(yc,[],3),[],9);
            D = @(i,j) By(:,(j-1)*3+i); % gets partial_i y^j
            cof = [
                D(2,2).*D(3,3)-D(3,2).*D(2,3),...%C(1,1)
               -D(1,2).*D(3,3)+D(3,2).*D(1,3),...%C(1,2)
                D(1,2).*D(2,3)-D(2,2).*D(1,3),...%C(1,3)
               -D(2,1).*D(3,3)+D(3,1).*D(2,3),...%C(2,1)
                D(1,1).*D(3,3)-D(3,1).*D(1,3),...%C(2,2)
               -D(1,1).*D(2,3)+D(2,1).*D(1,3),...%C(2,3)
                D(2,1).*D(3,2)-D(3,1).*D(2,2),...%C(3,1)
               -D(1,1).*D(3,2)+D(3,1).*D(1,2),...%C(3,2)
                D(1,1).*D(2,2)-D(2,1).*D(1,2)... %C(3,3)
                ];

            det = sum(By(:,1:3).*cof(:,1:3),2);
            if doDerivative

                D1 = Mesh.dx1; D2 = Mesh.dx2; D3 = Mesh.dx3;

                dDet = [sdiag(cof(:,1))*D1 + sdiag(cof(:,2))*D2 + sdiag(cof(:,3))*D3,...
                        sdiag(cof(:,4))*D1 + sdiag(cof(:,5))*D2 + sdiag(cof(:,6))*D3,...
                        sdiag(cof(:,7))*D1 + sdiag(cof(:,8))*D2 + sdiag(cof(:,9))*D3];
                
            end
        
        else
            
            dx1 = Mesh.mfdx1; dx2 = Mesh.mfdx2;  dx3 = Mesh.mfdx3;
            
            cof11 = @(yc)  dx2.D(yc(:,2)).*dx3.D(yc(:,3)) - dx3.D(yc(:,2)).*dx2.D(yc(:,3));
            cof12 = @(yc) -dx1.D(yc(:,2)).*dx3.D(yc(:,3)) + dx3.D(yc(:,2)).*dx1.D(yc(:,3));
            cof13 = @(yc)  dx1.D(yc(:,2)).*dx2.D(yc(:,3)) - dx2.D(yc(:,2)).*dx1.D(yc(:,3));
            cof21 = @(yc) -dx2.D(yc(:,1)).*dx3.D(yc(:,3)) + dx3.D(yc(:,1)).*dx2.D(yc(:,3));
            cof22 = @(yc)  dx1.D(yc(:,1)).*dx3.D(yc(:,3)) - dx3.D(yc(:,1)).*dx1.D(yc(:,3));
            cof23 = @(yc) -dx1.D(yc(:,1)).*dx2.D(yc(:,3)) + dx2.D(yc(:,1)).*dx1.D(yc(:,3));
            cof31 = @(yc)  dx2.D(yc(:,1)).*dx3.D(yc(:,2)) - dx3.D(yc(:,1)).*dx2.D(yc(:,2));
            cof32 = @(yc) -dx1.D(yc(:,1)).*dx3.D(yc(:,2)) + dx3.D(yc(:,1)).*dx1.D(yc(:,2));
            cof33 = @(yc)  dx1.D(yc(:,1)).*dx2.D(yc(:,2)) - dx2.D(yc(:,1)).*dx1.D(yc(:,2));
            
            det = dx1.D(yc(:,1)).*cof11(yc) + dx2.D(yc(:,1)).*cof12(yc) + dx3.D(yc(:,1)).*cof13(yc);
            if doDerivative

                dDet.dDet = @(x)  cof11(yc).*dx1.D(x(:,1)) + cof12(yc).*dx2.D(x(:,1)) + cof13(yc).*dx3.D(x(:,1)) +...
                                  cof21(yc).*dx1.D(x(:,2)) + cof22(yc).*dx2.D(x(:,2)) + cof23(yc).*dx3.D(x(:,2)) +...
                                  cof31(yc).*dx1.D(x(:,3)) + cof32(yc).*dx2.D(x(:,3)) + cof33(yc).*dx3.D(x(:,3));

                dDet.dDetadj = @(x) [ ...
                                 dx1.Dadj(cof11(yc).*x) + dx2.Dadj(cof12(yc).*x) + dx3.Dadj(cof13(yc).*x);...
                                 dx1.Dadj(cof21(yc).*x) + dx2.Dadj(cof22(yc).*x) + dx3.Dadj(cof23(yc).*x);...
                                 dx1.Dadj(cof31(yc).*x) + dx2.Dadj(cof32(yc).*x) + dx3.Dadj(cof33(yc).*x)
                                 ]; 
            
            
            end
            
        end
        
        
    otherwise
        error('only dim=2,3.');
end
function D = sdiag(v)
D =  spdiags(v(:),0,numel(v),numel(v));

function runMinimalExample
omega = [0 3 0 5]; m = [2 4]; type = 1;
Mesh  = TriMesh1(omega,m);
xc    = Mesh.xn + 1e-1*randn(size(Mesh.xn));

[det,dD]   = feval(mfilename,Mesh,xc);

fctn = @(x)getDeterminant(Mesh,x);
checkDerivative(fctn,xc(:));

fig = figure(1); clf;
subplot(1,2,1);
trisurf(Mesh.tri,xc(:,1),xc(:,2),0*xc(:,1),det); view([0 90]);
axis(omega);
title(sprintf('sum(V)=%f',sum(det)));

omega = [0 3 0 5 0 2]; m = [2 4 3];
Mesh  = TetraMesh1(omega,m);
xc    = Mesh.xn + 1e-1*randn(size(Mesh.xn));

[det,dD]   = getDeterminant(Mesh,xc);

fctn = @(x)getDeterminant(Mesh,x);
checkDerivative(fctn,xc(:));

figure(fig)
subplot(1,2,2);
tetramesh(Mesh.tri,xc,det,'FaceAlpha',.3);

axis(omega);
title(sprintf('sum(V)=%f',sum(det)));
