%==============================================================================
% This code is part of the VAMPIRE app for the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR/VAMPIRE.m 
%==============================================================================
% PCG solver for VAMPIRE
% 

function dy = FEMPIREsolveGN_PCG(rhs, H, maxIterCG, tolCG)

h = (H.omega(2:2:end)-H.omega(1:2:end))./H.m;
hd = prod(h);
% set pointers to simplify notation
Mesh    = H.Mesh;
P       = H.d2D.P;
Jac     = H.d2D.Jac;
dJac    = H.d2D.dJac;
dJacadj = H.d2D.dJacadj;
Tc      = H.d2D.Tc;
dT      = H.d2D.dT;
yc      = H.d2S.yc;
dres    = H.d2D.dres;
dim     = Mesh.dim;

% Tcmod(y) = T(P*y) * Jac(y)
% ==> dTcmod(y)  = Jac(y) * (dT(y)*P) + T(P*y) * dJac(y)
% ==> dTcmod'(y) = (P'* dT') Jac(y) + (dJac(y))' * Tc
dTcmod =    @(x) vecmatprod1(Jac,dT,H.d2D.P(reshape(x,[],dim))) + Tc .* dJac(reshape(x,[],dim));
dTcmodadj = @(x) H.d2D.P(vecmatprod1adj(dT,Jac,x)) + dJacadj(Tc.*x);

M         = @(x) dTcmodadj(dres'*H.d2D.d2psi*(dres * dTcmod(x)));
Hoperator = @(x) M(x) + H.d2S.d2S(x,H.omega,H.m,H.d2S.yc);

D         = H.d2D.diag(Mesh,yc)  +  H.d2S.diag(H.d2S.yc);
Preconditioner = @(x) D.\x; % Jacobi preconditioner
[dy,flag,relres,iter] = pcg(Hoperator,rhs,tolCG,maxIterCG,Preconditioner);
       
       
function p= vecmatprod1(Jac,dI,x)
% We aim to compute diag(Jac) * [dI] * x
dim  = size(dI,2);
n    = length(Jac);

x = reshape(x,[],dim);
p = zeros(n,1);
for d=1:dim
    p = p + dI(:,d) .* x(:,d);
end
p = p.*Jac;

function p = vecmatprod1adj(dI,Jac,x)
% We aim to compute dI' * [Jac] * x
dim  = size(dI,2);
n = length(Jac);

p = zeros(n,dim);
for d=1:dim
    p(:,d) = dI(:,d) .* Jac.*x;
end
