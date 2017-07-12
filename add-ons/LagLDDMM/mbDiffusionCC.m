%==============================================================================
% This code is part of the Matlab-based toolbox 
% LagLDDDM - A Lagrangian Gauss--Newton--Krylov Solver for Mass- and 
%                        Intensity-Preserving Diffeomorphic Image Registration
% 
% For details and license info see 
% - https://github.com/C4IR/FAIR.m/tree/master/add-ons/LagLDDMM
%
% function [Sc,dS,d2S] = diffusionST(vc,omega,m,varargin)
%
% Matrix-based spatio-temporal diffusion regularization energy for vc, where
% vc is cell-centered
%
% S(v) = 0.5 * \int_{\omega} alpha(1)*v(x)'*A*v(x)+ alpha(2)*v(x)'*B*v(x) dx,
%
% where A is the spatial gradient operator and B the time derivative operator.
%
% Input:
%
%   vc          instationary velocity field (cell-centered)
%   omega       spatial domain
%   m           number of discretization points
%   varargin    optional parameters (see below)
%
% Optional Input:
%
%   tspan       time span (default: [0 1])
%   nt          number of time points for velocity (default: computed from input)
%
%
% Output:
%
%   Sc          current value  (0.5 * hd * vc'* A *vc + 0.5*dt*vc'*B*vc)
%   dS          derivative     (hd * vc'*A )
%   d2S         Hessian        A
% ==================================================================================

function [Sc,dS,d2S] = mbDiffusionST(vc,omega,m,varargin)
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


persistent A omegaOld mOld alphaOld

if ~exist('mOld','var'),     mOld = [];     end;
if ~exist('omegaOld','var'), omegaOld = []; end;
if ~exist('alphaOld','var'), alphaOld = []; end;

alpha       = 1;
for k=1:2:length(varargin), % overwrites default parameter
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;


    build =  isempty(mOld) || isempty(omegaOld) ...
        || length(mOld) ~= length(m) || length(omegaOld) ~= length(omega) ...
        || any(mOld ~= m) || any(omegaOld~=omega)|| any(alphaOld~=alpha(1));
    if build,
        fprintf('%s - rebuilding op\n',mfilename);
        A = getDiffusionMatrix(omega,m,alpha(1));
        A = A'*A;
        mOld = m; omegaOld = omega; alphaOld = alpha(1); 
    end;
    dS  = vc'*A;
    Sc  = 0.5*dS*vc;
    d2S = A;

function B = getDiffusionMatrix(omega,m,alpha)

h   = (omega(2:2:end)-omega(1:2:end))./m;
hd  = prod(h);
B   = sqrt(alpha.*hd)*getSpaceGradientMatrix(omega,m);

function A = getSpaceGradientMatrix(omega,m)
dim = length(omega)/2;
h   = (omega(2:2:end)- omega(1:2:end)) ./m ;
I = @(i)  speye(m(i));

% setup regularizer for y
%
% S(y) =  INT |Dy(.,t)|^2 + INT |d/dt y(x,.)|^2 dx
switch dim
    case 2
        % build discrete derivative operators
        D1=spdiags(ones(m(1),1)*[-1 1],0:1,m(1)-1,m(1)); D1=D1/(h(1));
        D2=spdiags(ones(m(2),1)*[-1 1],0:1,m(2)-1,m(2)); D2=D2/(h(2));
        
        % spatial regularization
        A = [kron(I(2),D1); kron(D2,I(1))];
    case 3
        % build discrete derivative operators
        D1=spdiags(ones(m(1),1)*[-1 1],0:1,m(1)-1,m(1)); D1=D1/(h(1));
        D2=spdiags(ones(m(2),1)*[-1 1],0:1,m(2)-1,m(2)); D2=D2/(h(2));
        D3=spdiags(ones(m(3),1)*[-1 1],0:1,m(3)-1,m(3)); D3=D3/(h(3));
        % build gradient (scalar)
        A = [  ...
            kron(I(3),kron(I(2),D1)); ...
            kron(I(3),kron(D2,I(1))); ...
            kron(D3,kron(I(2),I(1)))...
            ];
        
end
A = kron(speye(dim), A);                  % nD gradient
