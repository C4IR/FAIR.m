%==============================================================================
% This code is part of the Matlab-based toolbox LagLDDDM - A Lagrangian Gauss--
% Newton--Krylov Solver for Mass- and Intensity-Preserving Diffeomorphic Image
% Registration
%
% For details and license info see
% - https://github.com/C4IR/LagLDDMM%
%==============================================================================
%
% function [Sc,dS,d2S] = mfDiffusionCC(vc,omega,m,varargin)
%
% computes diffusion regularization energy for vc, where vc cell-centered
%
% S(u) = 0.5 * \int_{\omega} u(x)'*A'*A*u(x) dx,
%
% where A is the gradient operator with Neuman boundary conditions.
%
% Input:
%
%   vc          displacement or stationary velocity field (cell-centered)
%   omega       spatial domain
%   m           number of discretization points
%   varargin    optional parameters (see below)
%
%
% Output:
%
%   Sc          current value  (alpha*0.5 * hd * vc'* A'*A *vc)
%   dS          derivative     (alpha*hd * vc'*A'*A )
%   d2S   		Hessian        alpha*hd*A'*A
%  if ~matrixFree,  d2S is sparse matrix; else, d2S is struct; endif
%
% ==================================================================================

function [Sc,dS,d2S] = diffusionCC(vc,omega,m,varargin)
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

persistent A omegaOld mOld alphaOld

if ~exist('mOld','var'),     mOld = [];     end;
if ~exist('omegaOld','var'), omegaOld = []; end;
if ~exist('alphaOld','var'), alphaOld = []; end;

alpha       = 1;
for k=1:2:length(varargin) % overwrites default parameter
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;


d2S.regularizer = regularizer;
d2S.alpha  = alpha(1);
d2S.B      = @(omega,m)     getDiffusionMatrix( omega,m,alpha(1));
d2S.d2S    = @(u,omega,m)   diffusionOperator(u,omega,m,alpha(1));
d2S.diag   = @(omega,m)     diag(omega,m,alpha(1));
d2S.solver = @spectralPrecondPCG;
d2S.res    = vc;
dS         = d2S.d2S(vc,omega,m)';
Sc         = .5*dS*vc;

function B = getDiffusionMatrix(omega,m,alpha)

h   = (omega(2:2:end)-omega(1:2:end))./m;
hd  = prod(h);
B   = sqrt(alpha.*hd)*getSpaceGradientMatrix(omega,m);



% get diagonal of d2S (interesting in matrix free mode)
function D = diag(omega,m,alpha)
dim   = numel(omega)/2;
one   = @(i,j) One(omega,m,i,j);
hd    = prod((omega(2:2:end)-omega(1:2:end))./m);

if dim == 2
    Dx   = [ one(1,1) + one(1,2);
        one(2,1) + one(2,2)];
else
    Dx   = [ ...
        one(1,1)+one(1,2)+one(1,3);
        one(2,2)+one(2,1)+one(2,3);
        one(3,3)+one(3,1)+one(3,2)];
end;
D  = hd*alpha(1) * Dx;


% helper for computation of diag(d2S)
function o = One(omega,m,i,j)
h = (omega(2:2:end)-omega(1:2:end))./m;
o = ones(m)/h(j)^2;
switch j
    case 1, o(2:end-1,:,:) = 2*o(2:end-1,:,:);
    case 2, o(:,2:end-1,:) = 2*o(:,2:end-1,:);
    case 3, o(:,:,2:end-1) = 2*o(:,:,2:end-1);
end;
o = o(:);

% matrix free implementation of diffusion operator
function Ay = diffusionOperator(vc,omega,m,alpha)
dim = numel(omega)/2;
h   = (omega(2:2:end)-omega(1:2:end))./m;
hd  = prod(h(1:dim));

switch dim
    case 2
        d1 = @(Y) (Y(2:end,:)-Y(1:end-1,:))/h(1);
        d2 = @(Y) (Y(:,2:end)-Y(:,1:end-1))/h(2);
        
        d1T = @(Y) reshape([-Y(1,:);Y(1:end-1,:)-Y(2:end,:);Y(end,:)],[],1)/h(1);
        d2T = @(Y) reshape([-Y(:,1),Y(:,1:end-1)-Y(:,2:end),Y(:,end)],[],1)/h(2);
        
        
        vc   = reshape(vc,[m dim]);
        Ay = hd*  alpha * ....
            [ (d1T(d1(vc(:,:,1))) + d2T(d2(vc(:,:,1)))); ...
            (d1T(d1(vc(:,:,2))) + d2T(d2(vc(:,:,2))))];
        Ay = Ay(:);
        
    case 3
        d1 = @(Y) (Y(2:end,:,:)-Y(1:end-1,:,:))/h(1);
        d2 = @(Y) (Y(:,2:end,:)-Y(:,1:end-1,:))/h(2);
        d3 = @(Y) (Y(:,:,2:end)-Y(:,:,1:end-1))/h(3);
        
        d1T = @(Y) reshape(d1t(Y),[],1)/h(1);
        d2T = @(Y) reshape(d2t(Y),[],1)/h(2);
        d3T = @(Y) reshape(d3t(Y),[],1)/h(3);
        
        
        vc  = reshape(vc,[m dim]);
        Ay  = hd*  alpha * ....
            [ ...
            (d1T(d1(vc(:,:,:,1))) + d2T(d2(vc(:,:,:,1))) + d3T(d3(vc(:,:,:,1)))); ...
            (d1T(d1(vc(:,:,:,2))) + d2T(d2(vc(:,:,:,2))) + d3T(d3(vc(:,:,:,2)))); ...
            (d1T(d1(vc(:,:,:,3))) + d2T(d2(vc(:,:,:,3))) + d3T(d3(vc(:,:,:,3)))); ...
            ];
        Ay = Ay(:);
    otherwise
        error('%s - dimension %d not supported.',mfilename,dim);
end


% partial derivative operator for x1 (mf)
function y = d1t(Y)
m = size(Y);
y = zeros(m+[1,0,0]);
y(1,:,:) = -Y(1,:,:);
y(2:end-1,:,:) = Y(1:end-1,:,:)-Y(2:end,:,:);
y(end,:,:) = Y(end,:,:);
% partial derivative operator for x2 (mf)
function y = d2t(Y)
m = size(Y);
y = zeros(m+[0,1,0]);
y(:,1,:) = -Y(:,1,:);
y(:,2:end-1,:) = Y(:,1:end-1,:)-Y(:,2:end,:);
y(:,end,:) = Y(:,end,:);
% partial derivative operator for x3 (mf)
function y = d3t(Y)
m = size(Y); if length(m) == 2, m = [m,1]; end;
y = zeros(m+[0,0,1]);
y(:,:,1) = -Y(:,:,1);
y(:,:,2:end-1) = Y(:,:,1:end-1)-Y(:,:,2:end);
y(:,:,end) = Y(:,:,end);

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
        A = [kron(I(3),kron(I(2),D1)); ...
            kron(I(3),kron(D2,I(1))); ...
            kron(D3,kron(I(2),I(1)))];
end
A = kron(speye(dim), A);                  % nD gradient

