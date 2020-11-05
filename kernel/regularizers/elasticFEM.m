%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration.
% For details see
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function [Sc,dS,d2S] = elastic(uc,omega,m,varargin)
%
% computes (linear) elastic regularization energy for uc = yc - yRef
%
% =============================================================================

function [Sc,dS,d2S] = elasticFEM(uc,omega,m,varargin)

if nargin == 0
  help(mfilename);
  runMinimalExample;
  Sc = 'endOfMinimalExample';
  return;
end

persistent A omegaOld mOld alphaOld muOld lambdaOld

matrixFree  = 0;
alpha       = 1;
mu          = 2;
lambda      = 1;

for k=1:2:length(varargin), % overwrites default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

if isempty(alpha)
  error('regularization parameter alpha is empty')
end;

dim = length(omega)/2;
h   = (omega(2:2:end)-omega(1:2:end))./m;
hd  = prod((omega(2:2:end)-omega(1:2:end))./m);

if not(matrixFree), % matrix-based
  
  build = isempty(mOld) || isempty(omegaOld) ...
    || isempty(alphaOld) || isempty(lambdaOld) || isempty(muOld) ...
    || length(mOld) ~= length(m) || length(omegaOld) ~= length(omega) ...
    || any(mOld ~= m) || any(omegaOld~=omega) || (alpha ~= alphaOld)...
    || (mu ~= muOld)|| (lambda ~= lambdaOld);
  if build,
    mOld      = m;
    omegaOld  = omega;
    alphaOld  = alpha;
    muOld     = mu;
    lambdaOld = lambda;
    
    A = alpha*getElasticFEMmatrix(omega,m,mu,lambda);
  end;
  dS  = uc'*A;
  Sc  = 0.5*dS*uc;
  d2S = A;
  
else % matrix-free
  
  d2S.regularizer = mfilename;
  d2S.alpha = alpha;
  d2S.operator = @(u,omega,m) elasticFEMOperator(u,omega,m,mu,lambda);
  d2S.matrix   = @(omega,m)   elasticFEMMatrix(omega,m,mu,lambda);
  d2S.diag     = @(omega,m)   getDiag(omega,m,mu,lambda,alpha);
  
  d2S.d2S   = @(u,omega,m) alpha*elasticFEMOperator(u,omega,m,mu,lambda);
  
  dS        = d2S.d2S(uc,omega,m)';
  Sc        = 0.5*dS*uc;
end

% get diagonal of d2S (interesting in matrix free mode)
function D = getDiag(omega,m,mu,lambda,alpha)
dim = length(omega)/2;
h   = (omega(2:2:end)-omega(1:2:end))./m;
a   = mu + lambda/2;
b   = mu/2;


if dim == 2,
  
  oneL = [1;2*ones(m(2)-1,1);1]*[1,2*ones(1 ,m(1)-1,1),1];
  D =  [
    (a/(3*h(1)^2))*oneL + (b/(3*h(2)^2))*oneL,...
    (b/(3*h(1)^2))*oneL + (a/(3*h(2)^2))*oneL
    ];
  
  
else
  D   = [ ...
    reshape( (2*mu+lambda)*one(1,1)+(mu+lambda)*(one(1,2)+one(1,3)),[],1)
    reshape( (2*mu+lambda)*one(2,2)+(mu+lambda)*(one(2,1)+one(2,3)),[],1)
    reshape( (2*mu+lambda)*one(3,3)+(mu+lambda)*(one(3,1)+one(3,2)),[],1) ];
end;
D = alpha*prod(h) *reshape(D,[],1);

% matrix free implementaion of elastic operator
function Au = elasticFEMOperator(uc,omega,m,mu,lambda)
dim = length(omega)/2;
h   = (omega(2:2:end)-omega(1:2:end))./m;
a   = mu + lambda/2;
b   = mu/2;
c   = 3*mu/4;
d   = 3*lambda/4;

switch dim,
  case 2,
    % function s = getEP(uc,omega,m,mu,lambda)
    
    
    uc = reshape(uc,[m+1,dim]);
    
    % GRAD U
    du = zeros(m(1),m(2),8);
    du(:,:,1) = squeeze( uc(2:end,  1:end-1,1) - uc(1:end-1,1:end-1,1))/h(1);
    du(:,:,2) = squeeze( uc(2:end,  2:end,  1) - uc(1:end-1,2:end,  1))/h(1);
    du(:,:,3) = squeeze( uc(1:end-1,2:end,  1) - uc(1:end-1,1:end-1,1))/h(2);
    du(:,:,4) = squeeze( uc(2:end,  2:end,  1) - uc(2:end,  1:end-1,1))/h(2);
    du(:,:,5) = squeeze( uc(2:end,  1:end-1,2) - uc(1:end-1,1:end-1,2))/h(1);
    du(:,:,6) = squeeze( uc(2:end,  2:end,  2) - uc(1:end-1,2:end,  2))/h(1);
    du(:,:,7) = squeeze( uc(1:end-1,2:end,  2) - uc(1:end-1,1:end-1,2))/h(2);
    du(:,:,8) = squeeze( uc(2:end,  2:end,  2) - uc(2:end,  1:end-1,2))/h(2);
    
    % weighting
    Tdu = zeros(m(1),m(2),8);
    Tdu(:,:,1) = (2*a)*du(:,:,1) +   a *du(:,:,2) +    d *du(:,:,7) +    d *du(:,:,8);
    Tdu(:,:,2) =    a *du(:,:,1) +(2*a)*du(:,:,2) +    d *du(:,:,7) +    d *du(:,:,8);
    Tdu(:,:,3) = (2*b)*du(:,:,3) +   b *du(:,:,4) +    c *du(:,:,5) +    c *du(:,:,6);
    Tdu(:,:,4) =    b *du(:,:,3) +(2*b)*du(:,:,4) +    c *du(:,:,5) +    c *du(:,:,6);
    Tdu(:,:,5) =    c *du(:,:,3) +   c *du(:,:,4) + (2*b)*du(:,:,5) +    b *du(:,:,6);
    Tdu(:,:,6) =    c *du(:,:,3) +   c *du(:,:,4) +    b *du(:,:,5) + (2*b)*du(:,:,6);
    Tdu(:,:,7) =    d *du(:,:,1) +   d *du(:,:,2) + (2*a)*du(:,:,7) +    a *du(:,:,8);
    Tdu(:,:,8) =    d *du(:,:,1) +   d *du(:,:,2) +    a *du(:,:,7) + (2*a)*du(:,:,8);
    
    % GRAD'
    d1T = @(Y) [-Y(1,:);Y(1:end-1,:)-Y(2:end,:);Y(end,:)]/h(1);
    d2T = @(Y) [-Y(:,1),Y(:,1:end-1)-Y(:,2:end),Y(:,end)]/h(2);
    
    Au = zeros([m+1,2]);
    Au( :     ,1:end-1,1) =                       d1T(Tdu(:,:,1));
    Au( :     ,2:end  ,1) = Au( :     ,2:end,1) + d1T(Tdu(:,:,2));
    Au(1:end-1, :     ,1) = Au(1:end-1, :   ,1) + d2T(Tdu(:,:,3));
    Au(2:end  , :     ,1) = Au(2:end  , :   ,1) + d2T(Tdu(:,:,4));
    Au( :     ,1:end-1,2) =                       d1T(Tdu(:,:,5));
    Au( :     ,2:end  ,2) = Au( :     ,2:end,2) + d1T(Tdu(:,:,6));
    Au(1:end-1, :     ,2) = Au(1:end-1, :   ,2) + d2T(Tdu(:,:,7));
    Au(2:end  , :     ,2) = Au(2:end  , :   ,2) + d2T(Tdu(:,:,8));
    
    Au = prod(h)/6 * reshape(Au,[],1);
    
    %{
    peye = @(j,k) spdiags(ones(m(j),1),k,m(j),m(j)+1);
    dx   = @(j)   spdiags(ones(m(j),1)*[-1,1],0:1,m(j),m(j)+1)/h(j);
    
    GRAD = spalloc(8*prod(m),2*prod(m+1),16*prod(m));
    GRAD = [
      kron(peye(2,0),dx(1))
      kron(peye(2,1),dx(1))
      kron(dx(2),peye(1,0))
      kron(dx(2),peye(1,1))
      ];
    GRAD = kron(speye(2,2),GRAD);

    W0 =  [
      2*a   a   0   0   0   0   d   d
      a 2*a   0   0   0   0   d   d
      0   0 2*b   b   c   c   0   0
      0   0   b 2*b   c   c   0   0
      0   0   c   c 2*b   b   0   0
      0   0   c   c   b 2*b   0   0
      d   d   0   0   0   0 2*a   a
      d   d   0   0   0   0   a 2*a ];
    
    W = kron(W0,speye(prod(m),prod(m)));
    A = getElasticFEMmatrix(omega,m,mu,lambda);
%     A = prod(h)/6*GRAD'*W*GRAD;
    test1 = norm(GRAD*uc(:)-du(:))
    if test1>1e-13,
      keyboard
    end
    test2 = norm(W*GRAD*uc(:)-Tdu(:))
    if test2>1e-13,
      keyboard
    end
    test3 = norm(A*uc(:)-Au(:))
    if test3>1e-13,
      keyboard
    end
    
    %}
    return;
  otherwise, error('nyi');
end;


function runMinimalExample
% --- 2D ---

alpha  = 100;
mu     = 2;
lambda = 1;


omega = [0 2 0 3];
m     = [13,14];
xc    = getNodalGrid(omega,m);

A = alpha*getElasticFEMmatrix(omega,m,mu,lambda);

n = size(A,1);
uv = @(j) ((1:n)' == j );

B = zeros(n,n);

for k=1:n,
  
  [S,dS] = feval(mfilename,uv(k),omega,m,'alpha',alpha,'mu',mu,'lambda',lambda,'matrixFree',1);
  B(:,k) = dS';
end;

figure(2); clf;
spy((A-B)>1e-10)
title('difference btw matrix and operaror version')



omega = [0 2 0 3];
m     = [17 19];
xc    = getNodalGrid(omega,m);
uc    = 1e-2*randn(size(xc));

[mbS, mbdS, mbd2S] = feval(mfilename,uc,omega,m,'alpha',1,'mu',1,'lambda',0,'matrixFree',0);
[mfS, mfdS, mfd2S] = feval(mfilename,uc,omega,m,'alpha',1,'mu',1,'lambda',0,'matrixFree',1);
fprintf('Elastic potential : %f (mb) | %f (mf) \n' , mbS,mfS);
RE = norm(mbdS-mfdS)/norm(mbdS);
fprintf('rel.error dS (mb vs. mf) : %e \n' , RE);
figure(1)
subplot(1,2,1);
spy(mbd2S)
title(sprintf('%s-spy(d2S)-(2D)',mfilename));




