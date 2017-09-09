clear all;
% verify getdeterminant matrix free implementation

%omega = [0 3 0 5]; m = [2 4];
%Mesh  = TriMesh1(omega,m);
omega = [0 3 0 5 0 2]; m = [2 4 3];
Mesh  = TetraMesh1(omega,m);
xc    = Mesh.xn;
% rotate x-y plane by 45 degree
yc(:,1) = (1/sqrt(2))*(xc(:,1) - xc(:,2));
yc(:,2) = (1/sqrt(2))*(xc(:,1) + xc(:,2));
yc(:,3) = xc(:,3);
xc = yc;

[det,dD]   = getDeterminant(Mesh,xc,'matrixFree',0);

[detMF,dDMF]   = getDeterminant(Mesh,xc,'matrixFree',1);

% error in determinant value
err_det = norm(det-detMF);

% error in dDet
x = rand(size(xc));
err_dD = norm(dD*x(:) - dDMF.dDet(x));

% error in adjoint of dDet 
z = rand(Mesh.ntri,1);
err_dDT = norm(dD'*z - dDMF.dDetadj(z));


