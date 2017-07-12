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
% Matrix-free spatio-temporal diffusion regularization energy for vc, where
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
%   Sc          current value   (0.5 * hd * vc'* A *vc + 0.5*dt*vc'*B*vc)
%   dS          derivative      (hd * vc'*A )
%   d2S         Hessian, struct A
% ==================================================================================

function [Sc,dS,d2S] = mfDiffusionST(vc,omega,m,varargin)
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
d2S.B      = @(omega,m)     getDiffusionMatrixST(omega,tspan,m,nt,alpha);
d2S.d2S    = @(u,omega,m)   diffusionOperatorST(u,omega,tspan,m,nt,alpha);
d2S.diag   = @(omega,m)     diag(omega,tspan,m,nt,alpha);
d2S.solver = @spectralPrecondPCG;
d2S.res    = vc;
dS         = d2S.d2S(vc,omega,m)';
Sc         = .5*dS*vc;

function B = getDiffusionMatrixST(omega,tspan,m,nt,alpha)

h   = (omega(2:2:end)-omega(1:2:end))./m;
hd  = prod(h);
% compute time-stepsize
dt = abs(tspan(2)-tspan(1))/nt;

% build gradient matrix for one transformation
Bx =  getSpaceGradientMatrix(omega,m);
% a = sqrt(alpha(1).*hd.*dt)*ones(nt+1,1);
a = sqrt(alpha(1).*hd.*dt*[1/2;ones(nt-1,1);1/2]);
% apply spatial regularization to all transformations and sum up in
% time
Bx =  kron(sdiag(a(:)),Bx);

% get time regularization matrix
b = sqrt(alpha(2)*hd.*dt);
Bt = b * getTimeDiffusionMatrix(nt,dt,prod(m),length(omega)/2);

% build B
B = [Bx;Bt];


function D = sdiag(v)
	D = diag(sparse(v(:)));

% get diagonal of d2S (interesting in matrix free mode)
function D = diag(omega,tspan,m,nt,alpha)
dim   = numel(omega)/2;
one   = @(i,j) One(omega,m,i,j);
hd    = prod((omega(2:2:end)-omega(1:2:end))./m);
% compute time stepsize
dt = abs(tspan(1)-tspan(2))/nt;

if dim == 2
    Dx = [ one(1,1) + one(1,2);
           one(2,1) + one(2,2) ];
    Dx = kron(hd*dt*ones(nt+1,1),Dx);

    Dt = [1/2;ones(nt-1,1);1/2]/dt^2;

    Dt = kron(hd*dt*Dt , ones(prod(m)*dim,1));

    D  = alpha(1)*Dx + alpha(2)*Dt;

else
    Dx   = [ ...
        one(1,1)+one(1,2)+one(1,3);
        one(2,2)+one(2,1)+one(2,3);
        one(3,3)+one(3,1)+one(3,2)];
    Dx   = kron(hd*dts(:),Dx);

    Dt          = zeros(mTime+1,1);
    Dt(1:end-1) = (1./dt(:));
    Dt(2:end  ) = Dt(2:end) + (1./dt(:));

    Dt = kron(hd*Dt, ones(prod(mSpace)*dim,1));

    D  = alpha(1)*Dx + alpha(2)*Dt;

end;

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

% matrix free implementation of spatio-temporal diffusion operator
function Ay = diffusionOperatorST(vc,omega,tspan,m,nt,alpha)
dim = numel(omega)/2;
h   = (omega(2:2:end)-omega(1:2:end))./m;
hd  = prod(h(1:dim));
dt  = abs(tspan(1)-tspan(2))/nt;
w   = dt*[1/2;ones(nt-1,1);1/2];

switch dim
    case 2
        d1 = @(Y) (Y(2:end,:)-Y(1:end-1,:))/h(1);
        d2 = @(Y) (Y(:,2:end)-Y(:,1:end-1))/h(2);

        d1T = @(Y) reshape([-Y(1,:);Y(1:end-1,:)-Y(2:end,:);Y(end,:)],[],1)/h(1);
        d2T = @(Y) reshape([-Y(:,1),Y(:,1:end-1)-Y(:,2:end),Y(:,end)],[],1)/h(2);


        vc   = reshape(vc,[],nt+1);
        Ay = zeros(prod(m)*dim,nt+1);
        % spatial diffusion
        for k=1:nt+1
            vct = reshape(vc(:,k),[m dim]);
            Ay(:,k) = w(k) * hd*  alpha(1) * ....
                [ (d1T(d1(vct(:,:,1))) + d2T(d2(vct(:,:,1)))); ...
                (d1T(d1(vct(:,:,2))) + d2T(d2(vct(:,:,2))))];
        end
        Ay = Ay(:);

    case 3
        d1 = @(Y) (Y(2:end,:,:)-Y(1:end-1,:,:))/h(1);
        d2 = @(Y) (Y(:,2:end,:)-Y(:,1:end-1,:))/h(2);
        d3 = @(Y) (Y(:,:,2:end)-Y(:,:,1:end-1))/h(3);

        d1T = @(Y) reshape(d1t(Y),[],1)/h(1);
        d2T = @(Y) reshape(d2t(Y),[],1)/h(2);
        d3T = @(Y) reshape(d3t(Y),[],1)/h(3);


        vc   = reshape(vc,[],nt+1);
        Ay = zeros(prod(m)*dim,nt+1);
        % spatial diffusion
        for k=1:nt+1
            vct = reshape(vc(:,k),[m dim]);
            Ay(:,k) = w(k) * hd*  alpha(1) * ....
                [ ...
                (d1T(d1(vct(:,:,:,1))) + d2T(d2(vct(:,:,:,1))) + d3T(d3(vct(:,:,:,1)))); ...
                (d1T(d1(vct(:,:,:,2))) + d2T(d2(vct(:,:,:,2))) + d3T(d3(vct(:,:,:,2)))); ...
                (d1T(d1(vct(:,:,:,3))) + d2T(d2(vct(:,:,:,3))) + d3T(d3(vct(:,:,:,3)))); ...
                ];
        end
        Ay = Ay(:);
    otherwise
        error('%s - dimension %d not supported.',mfilename,dim);
end
% time diffusion
if alpha(2)>0
    Dt  = @(Y) (Y(:,2:end)-Y(:,1:end-1))/dt;
    DtT = @(Y) ([-Y(:,1),Y(:,1:end-1)-Y(:,2:end),Y(:,end)])/dt;
    At  = DtT(Dt(vc));
    Ay  = Ay + (alpha(2)* hd*dt)*At(:);
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

function A = getTimeDiffusionMatrix(nt,dt,n,dim)
Dt=spdiags(ones(nt+1,1)*[-1/dt 1/dt],0:1,nt,nt+1);
A = kron(Dt, speye(n*dim)); % 2nd deriveative over time for two components

