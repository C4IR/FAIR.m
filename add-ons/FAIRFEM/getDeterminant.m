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

doDerivative = (nargout>1);
for k=1:2:length(varargin)     % overwrites default parameter
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

matrixFree = regularizer('get','matrixFree');
dim = Mesh.dim;
yc = reshape(yc,[],dim);

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
        % compute derivative
        GRAD = Mesh.GRAD;
        By = reshape(GRAD*reshape(yc,[],3),[],9);
        D = @(i,j) By(:,(i-1)*3+j);
        cof = [
            D(2,2).*D(3,3)-D(2,3).*D(3,2),...
            D(2,1).*D(3,3)-D(2,3).*D(3,2),...
            D(2,1).*D(3,2)-D(2,2).*D(3,1),...
            ];

        det = sum(By(:,1:3).*cof(:,1:3),2);
        if doDerivative,
            
            D1 = Mesh.dx1; D2 = Mesh.dx2; D3 = Mesh.dx3;
            Z = sparse(size(D1,1), size(D2,2));
            D = @(i,j) sdiag(By(:,(i-1)*3+j));
            
            dCof{1} = [  Z                 , D(3,3)*D2-D(3,2)*D3, D(2,2)*D3-D(2,3)*D2];
            dCof{2} = [  Z                 , D(3,3)*D1-D(3,2)*D3, D(2,1)*D3-D(2,3)*D2];
            dCof{3} = [  Z                 , D(3,2)*D1-D(3,1)*D2, D(2,1)*D2-D(2,2)*D1];
            
            dDet = [sdiag(cof(:,1))*D1 + sdiag(cof(:,2))*D2 + sdiag(cof(:,3))*D3,Z,Z]...
                + sdiag(By(:,1))*dCof{1} + sdiag(By(:,2))*dCof{2} + sdiag(By(:,3))*dCof{3};
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
