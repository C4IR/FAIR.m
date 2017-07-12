%==============================================================================
% This code is part of the VAMPIRE app for the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR/VAMPIRE.m 
%==============================================================================
% PCG solver for VAMPIRE
% 

function dy = VAMPIREsolveGN_PCG(rhs, H, maxIterCG, tolCG)

h = (H.omega(2:2:end)-H.omega(1:2:end))./H.m;
hd = prod(h);
% set pointers to simplify notation
P       = H.d2D.P;
Jac     = H.d2D.Jac;
dJac    = H.d2D.dJac;
dJacadj = H.d2D.dJacadj;
Tc      = H.d2D.Tc;
dT      = H.d2D.dT;
yc      = H.d2S.yc;
dres    = H.d2D.dres;
w2y     = H.d2D.w2y;
y2w     = H.d2D.y2w;

% Tcmod(y) = T(P*y) * Jac(y)
% ==> dTcmod(y)  = Jac(y) * (dT(y)*P) + T(P*y) * dJac(y)
% ==> dTcmod'(y) = (P'* dT') Jac(y) + (dJac(y))' * Tc
dTcmod =    @(x) vecmatprod1(Jac,dT,P(x)) + Tc .* dJac(yc,H.m,x);
dTcmodadj = @(x) P(vecmatprod1adj(dT,Jac,x)) + dJacadj(yc,H.m,Tc.*x);

% preconditioning
Dmf = geometryMexC(H.d2S.yc,H.m,'VAMPIREdiag',Jac,dT(:),Tc,h);
D =  Dmf + H.d2S.diag(H.d2S.yc);
PC = @(x) y2w(D).\x; % Jacobi preconditioner

% Afctn = @(x) hd * dTcmodadj(dres'*(dres * dTcmod(x))) + H.d2S.d2S(x,H.omega,H.m,yc);
operator = @(x) y2w(hd * dTcmodadj(dres'*(dres * dTcmod(x))) + H.d2S.d2S(x,H.omega,H.m,yc));
Afctn = @(x) operator(w2y(x));

[dy,flag,relres,iter] = pcg(Afctn,rhs,tolCG,maxIterCG,PC);

function p= vecmatprod1(Jac,dI,x)
% We aim to compute diag(Jac) * [dI] * x
%
% size(Jac) = [n, 1]
% size(dI) = [nVol * n , dim * n]
% size(x)  = [dim*n,1]
dim  = size(dI,2);
nVol = size(dI,3);
n    = length(Jac);

x = reshape(x,[],dim);
p = zeros(n,nVol);
for vol=1:nVol,
    for d=1:dim,
        p(:,vol) = p(:,vol) + dI(:,d,vol) .* x(:,d);
    end
    p(:,vol) = p(:,vol).*Jac;
end
p = p(:);

function p = vecmatprod1adj(dI,Jac,x)
% We aim to compute dI' * [Jac] * x
%
% size(Jac) = [n,1]
% size(dI)  = [nVol*n, dim*n]
% size(x)   = [nVol*n,1]
dim  = size(dI,2);
nVol = size(dI,3);
n = length(Jac);

x = reshape(x,[],nVol);
p = zeros(n,dim);
for vol=1:nVol,
    for d=1:dim,
        p(:,d) = p(:,d) + dI(:,d,vol) .* x(:,vol);
    end
end
for d=1:dim,
    p(:,d) = p(:,d).*Jac;
end
p = p(:);