%==============================================================================
% This code is part of the Finite Element Method app for the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR/FAIRFEM 
%==============================================================================
%
% function [Dx1,Dx2,Dx3] = getGradientMatrixFEM(Mesh)
%
% builds gradient operator for tetrahedral finite element discretization
% with piecewise linear basis functions
%
% Input:
%
% Mesh          - description of mesh, struct
%
% Output:
%
% [Dx1,Dx2,Dx3] - discrete partial derivative operators
%
%==============================================================================


function [Dx1,Dx2,Dx3] = getGradientMatrixFEM(Mesh,matrixFree)
if nargin==0, help(mfilename); runMinimalExample; return; end;

if not(exist('matrixFree','var'))||isempty(matrixFree),
    matrixFree=0;
end

Dx3 = [];

xn  = Mesh.xn;
vol = Mesh.vol;
dim = Mesh.dim;

flag = [num2str(dim) 'D'];
if matrixFree
    flag = [flag '-mf'];
else
    flag = [flag '-mb'];
end

switch flag
    case '2D-mb'
        % compute edges
        e1 =  Mesh.P1*xn - Mesh.P3*xn;
        e2 =  Mesh.P2*xn - Mesh.P3*xn;
        % compute gradients of basis functions
        dphi1 =   [ e2(:,2) -e2(:,1)]./[2*vol,2*vol];
        dphi2 =   [-e1(:,2)  e1(:,1)]./[2*vol,2*vol];
        dphi3 = -dphi1-dphi2;
 
        Dx1 = sdiag(dphi1(:,1))*Mesh.P1 +  sdiag(dphi2(:,1))*Mesh.P2  + sdiag(dphi3(:,1))*Mesh.P3;
        Dx2 = sdiag(dphi1(:,2))*Mesh.P1 +  sdiag(dphi2(:,2))*Mesh.P2  + sdiag(dphi3(:,2))*Mesh.P3;
        
        if nargout==1, % return gradient operator for vector field
            B = [Dx1;Dx2];
            Dx1 = blkdiag(B,B);
        end
        
    case '2D-mf'
        % compute edges
        x1 =  Mesh.mfPi(xn,1);
        x2 =  Mesh.mfPi(xn,2);
        x3 =  Mesh.mfPi(xn,3);
        e1 =  x1 - x3;
        e2 =  x2 - x3;
        
        % compute gradients of basis functions
        dphi1 =   [ e2(:,2) -e2(:,1)]./[2*vol,2*vol];
        dphi2 =   [-e1(:,2)  e1(:,1)]./[2*vol,2*vol];
        dphi3 = -dphi1 - dphi2;
        
        Dxi    = @(x,i)   dphi1(:,i).*Mesh.mfPi(x,1) ...
                        + dphi2(:,i).*Mesh.mfPi(x,2) ...
                        + dphi3(:,i).*Mesh.mfPi(x,3);
                    
        Dxiadj = @(x,i)   Mesh.mfPi(dphi1(:,i).*x,1) ....
                        + Mesh.mfPi(dphi2(:,i).*x,2) ...
                        + Mesh.mfPi(dphi3(:,i).*x,3);
        if nargout==1
            % matrix free mode, return D*xc and operators
            Dx1.D    = @(x) getGradientMF(x,Dxi,dim,'Du');
            Dx1.Dadj = @(x) getGradientMF(x,Dxiadj,dim,'DTu');
        else
            Dx1.D    = @(x) Dxi(x,1);
            Dx1.Dadj = @(x) Dxiadj(x,1);
            Dx2.D    = @(x) Dxi(x,2);
            Dx2.Dadj = @(x) Dxiadj(x,2);
        end
    case '3D-mf'
        % compute edges
        v1   = Mesh.mfPi(xn,1);
        v2   = Mesh.mfPi(xn,2);
        v3   = Mesh.mfPi(xn,3);
        v4   = Mesh.mfPi(xn,4);
        e1   =  v1 - v4;
        e2   =  v2 - v4;
        e3   =  v3 - v4;
        % compute inverse transformation to reference element
        cofA =  [
              e2(:,2).*e3(:,3)-e2(:,3).*e3(:,2), ...
            -(e1(:,2).*e3(:,3)-e1(:,3).*e3(:,2)), ...
              e1(:,2).*e2(:,3)-e1(:,3).*e2(:,2), ...
            -(e2(:,1).*e3(:,3)-e2(:,3).*e3(:,1)), ...
              e1(:,1).*e3(:,3)-e1(:,3).*e3(:,1), ...
            -(e1(:,1).*e2(:,3)-e1(:,3).*e2(:,1)), ...
              e2(:,1).*e3(:,2)-e2(:,2).*e3(:,1), ...
            -(e1(:,1).*e3(:,2)-e1(:,2).*e3(:,1)), ...
              e1(:,1).*e2(:,2)-e1(:,2).*e2(:,1), ...
            ];
        detA = e1(:,1).*cofA(:,1) + e2(:,1).*cofA(:,2)+ e3(:,1).*cofA(:,3);
        
        
        % compute gradients of basis functions
         dphi1 =   cofA(:,[1 4 7])./repmat(detA,[1 3]);
        dphi2 =   cofA(:,[2 5 8])./repmat(detA,[1 3]);
        dphi3 =   cofA(:,[3 6 9])./repmat(detA,[1 3]);
       dphi4 =   -dphi1 - dphi2 - dphi3;
        
        Dxi    = @(x,i) dphi1(:,i).*Mesh.mfPi(x,1) + dphi2(:,i).*Mesh.mfPi(x,2) ...
            + dphi3(:,i).*Mesh.mfPi(x,3) + dphi4(:,i).*Mesh.mfPi(x,4);
        Dxiadj = @(x,i) Mesh.mfPi(dphi1(:,i).*x,1) + Mesh.mfPi(dphi2(:,i).*x,2) ...
            + Mesh.mfPi(dphi3(:,i).*x,3) + Mesh.mfPi(dphi4(:,i).*x,4);
        
        if nargout==1
            % matrix free mode, return D*xc and operators
            Dx1.D    = @(x) getGradientMF(x,Dxi,dim,'Du');
            Dx1.Dadj = @(x) getGradientMF(x,Dxiadj,dim,'DTu');
        else
            Dx1.D    = @(x) Dxi(x,1);
            Dx1.Dadj = @(x) Dxiadj(x,1);
            Dx2.D    = @(x) Dxi(x,2);
            Dx2.Dadj = @(x) Dxiadj(x,2);
            Dx3.D    = @(x) Dxi(x,3);
            Dx3.Dadj = @(x) Dxiadj(x,3);
        end
        
        
        
    case '3D-mb'
        % compute edges
        e1   =  Mesh.P1*xn - Mesh.P4*xn;
        e2   =  Mesh.P2*xn - Mesh.P4*xn;
        e3   =  Mesh.P3*xn - Mesh.P4*xn;
        % compute inverse transformation to reference element
        cofA =  [
            e2(:,2).*e3(:,3)-e2(:,3).*e3(:,2), ...
            -(e1(:,2).*e3(:,3)-e1(:,3).*e3(:,2)), ...
            e1(:,2).*e2(:,3)-e1(:,3).*e2(:,2), ...
            -(e2(:,1).*e3(:,3)-e2(:,3).*e3(:,1)), ...
            e1(:,1).*e3(:,3)-e1(:,3).*e3(:,1), ...
            -(e1(:,1).*e2(:,3)-e1(:,3).*e2(:,1)), ...
            e2(:,1).*e3(:,2)-e2(:,2).*e3(:,1), ...
            -(e1(:,1).*e3(:,2)-e1(:,2).*e3(:,1)), ...
            e1(:,1).*e2(:,2)-e1(:,2).*e2(:,1), ...
            ];
        detA = e1(:,1).*cofA(:,1) + e2(:,1).*cofA(:,2)+ e3(:,1).*cofA(:,3);
        
        
        % compute gradients of basis functions
        dphi1 =   cofA(:,[1 4 7])./repmat(detA,[1 3]);
        dphi2 =   cofA(:,[2 5 8])./repmat(detA,[1 3]);
        dphi3 =   cofA(:,[3 6 9])./repmat(detA,[1 3]);
        dphi4 =   - dphi1 - dphi2 - dphi3;
        
        Dxi = @(i) sdiag(dphi1(:,i))*Mesh.P1 +  sdiag(dphi2(:,i))*Mesh.P2  ...
            + sdiag(dphi3(:,i))*Mesh.P3 +  sdiag(dphi4(:,i))*Mesh.P4;
        
        Dx1 = Dxi(1); Dx2 = Dxi(2); Dx3 = Dxi(3);
        
        if nargout==1, % return gradient operator for vector field
            B = [Dx1;Dx2;Dx3];
            Dx1 = blkdiag(B,B,B);
        end
        
    otherwise
        error('unknown flag: %s',flag);
end

function Du = getGradientMF(x,Dxi,dim,flag)
flag = [num2str(dim) 'D-' flag];
switch flag
    case '2D-Du'        
        x = reshape(x,[],dim);
        Du  = [  Dxi(x(:,1),1);  Dxi(x(:,1),2);  Dxi(x(:,2),1);  Dxi(x(:,2),2); ];
    case '2D-DTu'
        x = reshape(x,[],4);
        Du  =  [Dxi(x(:,1),1)+Dxi(x(:,2),2); Dxi(x(:,3),1)+Dxi(x(:,4),2)];
    case '3D-Du'        
        x = reshape(x,[],dim);
        Du  = [  
                Dxi(x(:,1),1);  Dxi(x(:,1),2);  Dxi(x(:,1),3);...
                Dxi(x(:,2),1);  Dxi(x(:,2),2);  Dxi(x(:,2),3);...
                Dxi(x(:,3),1);  Dxi(x(:,3),2);  Dxi(x(:,3),3);...
               ];
    case '3D-DTu'
        x = reshape(x,[],9);
        Du  =  [ 
                    Dxi(x(:,1),1)+Dxi(x(:,2),2)+Dxi(x(:,3),3);...
                    Dxi(x(:,4),1)+Dxi(x(:,5),2)+Dxi(x(:,6),3);...
                    Dxi(x(:,7),1)+Dxi(x(:,8),2)+Dxi(x(:,9),3);...
                ];
    otherwise
        error('nyi');
end


% shortcut for sparse diagonal matrices
function D = sdiag(v)
D = spdiags(v(:),0,numel(v),numel(v));

function runMinimalExample
omega = [0 1 0 3]; m = [ 3 5];
Mesh = TriMesh1(omega,m);
xn = Mesh.xn + 6*1e-2 * randn(size(Mesh.xn));

D = feval(mfilename,Mesh,0);

figure(1); clf;
subplot(1,3,1);
triplot(Mesh.tri,xn(:,1),xn(:,2)); axis(omega); hold on;
plot(xn(:,1),xn(:,2),'or');
title('discretization');
subplot(1,3,2); spy(D);
title('gradient operator B');
subplot(1,3,3); spy(D'*D);
title('Laplace operator B''*B');


Dmf = feval(mfilename,Mesh,1);

err1 = norm(D*xn(:)-Dmf.D(xn))/norm(D*xn(:));
x = randn(size(xn)); y=randn(4*size(Mesh.tri,1),1);
err2 = abs(x(:)'*Dmf.Dadj(y)-y(:)'*Dmf.D(x));
fprintf('MFerror: %1.3e;  AdjointError: %1.3e OK? %d\n',err1,err2,max(err1,err2)<1e-12);
%========= 3D ==========
omega = [0 1 0 3 0 3]; m = [ 3 5 4];
Mesh = getTriangleMesh(omega,m,'matrixFree',0);
xn = Mesh.xn + 6*1e-2 * randn(size(Mesh.xn));

D = feval(mfilename,Mesh);

figure(1); clf;
subplot(1,3,1);
tetramesh(Mesh.tri,Mesh.xn,'FaceAlpha',.3); hold on;
plot3(xn(:,1),xn(:,2),xn(:,3),'or');
title('discretization');
subplot(1,3,2); spy(D);
title('gradient operator B');
subplot(1,3,3); spy(D'*D);
title('Laplace operator B''*B');

MeshMF = getTriangleMesh(omega,m,'matrixFree',1);
Dmf = feval(mfilename,MeshMF);

err1 = norm(D*xn(:)-Dmf.D(xn))/norm(D*xn(:));
x = randn(size(xn)); y=randn(9*size(Mesh.tri,1),1);
err2 = abs(x(:)'*Dmf.Dadj(y)-y(:)'*Dmf.D(x));
fprintf('MFerror: %1.3e AdjointError: %1.3e OK? %d\n',err1,err2,max(err1,err2)<1e-12);






