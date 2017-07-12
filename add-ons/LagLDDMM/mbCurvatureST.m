%==============================================================================
% This code is part of the Matlab-based toolbox 
% LagLDDDM - A Lagrangian Gauss--Newton--Krylov Solver for Mass- and 
%                        Intensity-Preserving Diffeomorphic Image Registration
% 
% For details and license info see 
% - https://github.com/C4IR/FAIR.m/tree/master/add-ons/LagLDDMM
%
% function [Sc,dS,d2S] = curvatureST(uc,omega,m,varargin)
%
% Matrix-free spatio-temporal curvature regularization energy for vc
% where vc is cell-centered
%
% Sv) = 0.5 * \int_{\omega}\int_0^1 
%               alpha(1)*v(x,t)'*A*v(x,t)+ alpha(2)*v(x,t)'*B*v(x,t) dx dt,
%
% where A is the curvature operator and B the first-order time derivative
% operator.
%
% Input:
%
%   vc          instationary velocity field (cell-centered)
%   omega       spatial domain
%   m           number of discretization points in space
%   varargin    optional parameters (see below)
%
% Optional Input:
%   tspan       time span (default: [0 1])
%   nt          number of time points for velocity
%
%
% Output:
%
%   Sc          current value  (0.5 * hd * uc'* A *uc)
%   dS          derivative     (hd * uc'*A )
%   d2S         Hessian        A
%  if ~matrixFree,  d2S is sparse matrix; else, d2S is struct; endif
%
% see also curvature.m
% ==================================================================================

function [Sc,dS,d2S] = mbCurvatureST(vc,omega,m,varargin)
if nargin == 0
    help(mfilename);
    return;
end

if strcmp(vc,'para')
    Sc = 'cell-centered';       % grid
    dS = 0;                     % matrixFree
    d2S = 'pcg';  % solver
    return;
end


persistent A omegaOld mOld alphaOld ntOld tspanOld

if ~exist('mOld','var'),     mOld = [];     end;
if ~exist('omegaOld','var'), omegaOld = []; end;
if ~exist('alphaOld','var'), alphaOld = []; end;

alpha       = [1 1e-3];
tspan       = [0 1];
nt          = [];
for k=1:2:length(varargin), % overwrites default parameter
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

dim = numel(omega)/2;

if isempty(nt), % roughly estimate nt
    nt = round(numel(vc)/(prod(m)*dim))-1;
end


    build =  isempty(mOld) || isempty(omegaOld) ...
        || length(mOld) ~= length(m) || length(omegaOld) ~= length(omega) ...
        || any(mOld ~= m) || any(omegaOld~=omega)|| any(alphaOld~=alpha) || any(tspanOld~=tspan) || any(ntOld~=nt);
    if build,
        fprintf('%s - rebuilding op\n',mfilename);
        A = getCurvatureMatrixST(omega,tspan,m,nt,alpha);
        A = A'*A;
        mOld = m; omegaOld = omega; alphaOld = alpha; ntOld = nt; tspanOld = tspan;
    end;
    dS  = vc'*A;
    Sc  = 0.5*dS*vc;
    d2S = A;


function C = getCurvatureMatrixST(omega,tspan,m,nt,alpha)

h   = (omega(2:2:end)-omega(1:2:end))./m;
hd  = prod(h);
% compute time-stepsize
dt = abs(tspan(2)-tspan(1))/nt;

% build gradient matrix for one transformation
C =  getCurvatureMatrix(omega,m);
a = sqrt(alpha(1).*hd.*dt*[1/2;ones(nt-1,1);1/2]);
% apply spatial regularization to all transformations and sum up in
% time
C =  kron(sdiag(a(:)),C);

% get time regularization matrix
if alpha(2)>0
    b = sqrt(alpha(2)*hd.*dt);
    Bt = b * getTimeDiffusionMatrix(nt,dt,prod(m),length(omega)/2);
    C  = [C;Bt];
end

function D = sdiag(v)
D = diag(sparse(v(:)));

function A = getTimeDiffusionMatrix(nt,dt,n,dim)
Dt=spdiags(ones(nt+1,1)*[-1/dt 1/dt],0:1,nt,nt+1);
A = kron(Dt, speye(n*dim)); % 2nd deriveative over time for two components
