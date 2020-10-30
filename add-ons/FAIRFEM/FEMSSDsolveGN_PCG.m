%==============================================================================
% This code is part of the VAMPIRE app for the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR/VAMPIRE.m 
%==============================================================================
% PCG solver for FEMSSD
% 

function dy = FEMSSDsolveGN_PCG(rhs, H, maxIterCG, tolCG)

% set pointers to simplify notation
Mesh    = H.Mesh;
P       = H.d2D.P;
Tc      = H.d2D.Tc;
dT      = H.d2D.dT;
yc      = H.d2S.yc;
dres    = H.d2D.dres;
dim     = Mesh.dim;

% Tcmod(y) = T(P*y)
% ==> dTcmod(y)  = (dT(y)*P)
% ==> dTcmod'(y) = (P'* dT')(y)
%dTcmod =    @(x) vecmatprod1(dT,H.d2D.P(reshape(x,[],dim)));
%dTcmodadj = @(x) H.d2D.P(vecmatprod1adj(dT,x));

dTcmod =    @(x) vecmatprod1(dT,H.d2D.P(reshape(x,[],dim)));
dTcmodadj = @(x) H.d2D.P(dT.*x);

M         = @(x) dTcmodadj(dres'*H.d2D.d2psi*(dres * dTcmod(x)));
Hoperator = @(x) M(x) + H.d2S.d2S(x);

D         = H.d2D.diag()  +  H.d2S.diag(H.d2S.yc);
Preconditioner = @(x) D.\x; % Jacobi preconditioner
[dy,flag,relres,iter] = pcg(Hoperator,rhs,tolCG,maxIterCG,Preconditioner);
       
       
function p= vecmatprod1(dI,x)
% We aim to compute [dI] * x
dim  = size(dI,2);

x = reshape(x,[],dim);
p = sum(dI.*x,2);

