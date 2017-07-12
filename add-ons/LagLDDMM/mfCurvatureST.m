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

function [Sc,dS,d2S] = mfCurvatureST(vc,omega,m,varargin)
if nargin == 0
    help(mfilename);
    return;
end

if strcmp(vc,'para')
    Sc = 'cell-centered';       % grid
    dS = 1;                     % matrixFree
    d2S = @spectralPrecondPCG;  % solver
    return;
end


persistent A omegaOld mOld alphaOld ntOld tspanOld

if ~exist('mOld','var'),     mOld = [];     end;
if ~exist('omegaOld','var'), omegaOld = []; end;
if ~exist('alphaOld','var'), alphaOld = []; end;

matrixFree  = 0;
alpha       = [1 1e-3];
tspan       = [0 1];
nt          = [];
for k=1:2:length(varargin) % overwrites default parameter
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

dim = numel(omega)/2;

if isempty(nt) % roughly estimate nt
    nt = round(numel(vc)/(prod(m)*dim))-1;
end


d2S.regularizer = regularizer;
d2S.alpha  = alpha;
d2S.B      = @(omega,m)     getCurvatureMatrixST( omega,tspan,m,nt,alpha);
d2S.d2S    = @(u,omega,m)   curvatureOperatorST(u,omega,tspan,m,nt,alpha);
d2S.diag   = @(omega,m)     getCurvatureDiag(omega,tspan,m,nt,alpha); %getDiag(omega,tspan,m,nt,alpha);
d2S.solver = @spectralPrecondPCG;
d2S.res    = vc;
dS         = d2S.d2S(vc,omega,m)';
Sc         = .5*dS*vc;

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

% matrix free implementation of spatio-temporal curvature operator
function Ay = curvatureOperatorST(uc,omega,tspan,m,nt,alpha)
dim = numel(omega)/2;
h   = (omega(2:2:end)-omega(1:2:end))./m;
hd  = prod(h(1:dim));
dt  = abs(tspan(1)-tspan(2))/nt;
w   = dt*[1/2;ones(nt-1,1);1/2];
n   = prod(m);

uc   = reshape(uc,dim*n,[]);
ntp1 = size(uc,2);
Ay = 0*uc;

switch dim
    case 2
        D2x =  D2(1,omega,m);
        D2y =  D2(2,omega,m);
        
        % space: curvature
        for k=1:ntp1
            uct = reshape(uc(:,k),[m dim]);
            t1  = D2x*uct(:,:,1) + uct(:,:,1)*D2y;
            t1  = D2x*t1 + t1*D2y;
            t2  = D2x*uct(:,:,2) + uct(:,:,2)*D2y;
            t2  = D2x*t2 + t2*D2y;
            
            Ay(:,k) = (w(k)*hd*alpha(1))*[t1(:);t2(:)];
        end
        Ay = Ay(:);
    case 3
        % space: curvature
        d2  = @(i) D2(i,omega,m); % this is a shortcut to the discrete second derivative
        
        % the following line is a shortcut for
        %  - permuting the 3d-array using the permutation J
        %  - reshape it to a 2D-array of size q-by-prod(m)/q, where q=m(J(1))
        %  - multiply by A (which is q-by-q)
        %  - undo the reshape, i.e. make the result to m(J(1))-by-m(J(2))-by-m(J(3))
        %  - undo the permutation
        operate = @(A,z,J) ipermute(reshape(A*reshape(permute(z,J),m(J(1)),[]),m(J)),J);
        for k=1:ntp1 % loop over all time steps
            uct = reshape(uc(:,k),[m,dim]);
            for ell=1:3 % run over all components y^ell of Y=(Y^1,Y^2,Y^3)
                % compute
                % (I_3\otimes I_2\otimes d2(1) + I_3\otimes d2(2)\otimes I_1 ...
                % + d2(3)\otimes I_2\otimes I_1) y^ell
                %                 utt = 0;
                for rep=1:2 % go twice to have the bi-Laplacian
                    for j=1:3
                        z = operate(d2(j),uct(:,:,:,ell),[j,setdiff(1:3,j)]);
                        Ay((ell-1)*n+(1:n),k) = Ay((ell-1)*n+(1:n),k) + reshape(z,[],1);
                    end;
                    if rep==1 % overwrite uct and reset Ay. Compare to t1,t2 in dim=2 implementation
                        uct(:,:,:,ell) = reshape(Ay((ell-1)*n+(1:n),k),m);
                        Ay((ell-1)*n+(1:n),k) = 0;
                    end
                end;
            end
            Ay(:,k) = w(k)*hd*alpha(1)*Ay(:,k);
        end
        Ay = Ay(:);
    otherwise
        error('%s - dimension %d not supported.',mfilename,dim);
end
% time: diffusion
if alpha(2)>0
    d3  = @(Y) (Y(:,2:end)-Y(:,1:end-1))/dt;
    d3T = @(Y) ([-Y(:,1),Y(:,1:end-1)-Y(:,2:end),Y(:,end)])/dt;
    At = d3T(d3(uc));
    Ay = Ay(:) + (alpha(2)* hd*dt)*At(:);
end

function D = getCurvatureDiag(omega,tspan,m,nt,alpha)
% just for testing - not very efficient
dim = numel(omega)/2;
h   = (omega(2:2:end)-omega(1:2:end))./m;
hd  = prod(h(1:dim));
dt  = abs(tspan(1)-tspan(2))/nt;

C = getCurvatureMatrix(omega,m);
D = full(diag(C'*C'));


D = hd*dt*alpha(1)*repmat(D(:),1,1);
if alpha(2)>0
    errory('nyi')
end


function D = D2(i,omega,m)
h = (omega(2:2:end)-omega(1:2:end))./m;
D = spdiags(ones(m(i),1)*[1,-2,1],-1:1,m(i),m(i))/h(i)^2;
D([1,end]) = -D([2,end-1]);


function A = getTimeDiffusionMatrix(nt,dt,n,dim)
Dt=spdiags(ones(nt+1,1)*[-1/dt 1/dt],0:1,nt,nt+1);
A = kron(Dt, speye(n*dim)); % 2nd deriveative over time for two components

