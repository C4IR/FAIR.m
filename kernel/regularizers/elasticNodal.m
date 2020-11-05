%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function [Sc,dS,d2S] = elasticNodal(uc,omega,m,varargin)
%
% % computes (linear) elastic regularization energy for uc = yc - yRef
% % (nodal)
% %
% % S(u) = 0.5 * \int_{\omega} u(x) ' * B' * B * u(x) dx,
% %
% % where B is the elastic operator, see getElasticMatrixNodal.m
% %
% % Input:
% % ------
% %   uc           displacement field (staggered)
% %   omega        spatial domain
% %   m            number of discretization points
% %   varargin     optional parameters (see below)
% %
% % Output:
% % -------
% %   Sc          current value  (0.5 * hd * uc'*B'*B*uc)
% %   dS          derivative     (hd * uc'*B'*B)
% %   d2S           Hessian        (B'*B)
% %  if ~matrixFree,  d2S is sparse matrix; else, d2S is struct; endif        
%
% =============================================================================
function [Sc,dS,d2S] = elasticNodal(uc,omega,m,varargin)

if nargin == 0
    help(mfilename);
    runMinimalExample;
    Sc = 'endOfMinimalExample';
    return;
end

persistent A omegaOld mOld alphaOld muOld lambdaOld

% if ~exist('mOld','var'),     mOld = [];     end;
% if ~exist('omegaOld','var'), omegaOld = []; end;

matrixFree  = 0;
alpha       = 1;
mu          = 1;
lambda      = 0;

for k=1:2:length(varargin), % overwrites default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

if isempty(alpha)
  error('regularization parameter alpha is empty')
end;

dim = length(omega)/2; 
h   = (omega(2:2:end)-omega(1:2:end))./m;
h = 1+0*h


hd  = prod((omega(2:2:end)-omega(1:2:end))./m);


uc = reshape(uc,[m+1,dim]);

d110 = squeeze( uc(2:end,  1:end-1,1) - uc(1:end-1,1:end-1,1))/h(1); 
d111 = squeeze( uc(2:end,  2:end,  1) - uc(1:end-1,2:end,  1))/h(1); 
d120 = squeeze( uc(1:end-1,2:end,  1) - uc(1:end-1,1:end-1,1))/h(2); 
d121 = squeeze( uc(2:end,  2:end,  1) - uc(2:end,  1:end-1,1))/h(2); 
d210 = squeeze( uc(2:end,  1:end-1,2) - uc(1:end-1,1:end-1,2))/h(1); 
d211 = squeeze( uc(2:end,  2:end,  2) - uc(1:end-1,2:end,  2))/h(1); 
d220 = squeeze( uc(1:end-1,2:end,  2) - uc(1:end-1,1:end-1,2))/h(2); 
d221 = squeeze( uc(2:end,  2:end,  2) - uc(2:end,  1:end-1,2))/h(2); 

E = 4/3*(d110.^2 + d110.*d111 + d111.^2) ...
	+ 2/3*(d210.^2 + d210.*d211 + d211.^2) ...
	+ 2/3*(d120.^2 + d120.*d121 + d121.^2) ...
	+ 4/3*(d220.^2 + d220.*d221 + d221.^2) ...
	+ d120.*d210 + d120.*d211 +d210.*d121 + d211.*d121;
  
E0 = sum(E(:))*hd



A0 = [
  4/3  2/3 0   0   0   0   0   0
  2/3  4/3 0   0   0   0   0   0     
  0    0   2/3 1/3 1/2 1/2 0   0		
  0    0   1/3 2/3 1/2 1/2 0   0		
  0    0   1/2 1/2 2/3 1/3 0   0 		
  0    0   1/2 1/2 1/3 2/3 0   0 		
  0    0   0   0   0   0   4/3 2/3 		
  0    0   0   0   0   0   2/3 4/3  		
];

bigA = kron(A0,speye(prod(m),prod(m)));

G = [d110(:);d111(:);d120(:);d121(:);d210(:);d211(:);d220(:);d221(:)];

E1 = G'*bigA*G;

I  = @(j,k) spdiags(ones(m(j),1),k,m(j),m(j)+1);
dx = @(j) spdiags(ones(m(j),1)*[-1,1],0:1,m(j),m(j)+1) / h(j);

grad = kron(speye(2,2),[ 
  kron(I(2,0),dx(1))
  kron(I(2,1),dx(1))
  kron(dx(2),I(1,0))
  kron(dx(2),I(1,1))
]);


A = grad'*bigA*grad;

A = A(1:prod(m+1),1:prod(m+1));
figure(1); spy(A)

a2 = @(j) 4*ones(m(j)+1,1) - 2*((1:m(j)+1)' == 1 | (1:m(j)+1)' == m(j)+1 ); 
T = kron(spdiags(a2(2),0,m(2)+1,m(2)+1),spdiags(ones(m(1)+1,1)*[-1 2 -1],-1:1,m(1)+1,m(1)+1));

figure(2), 
spy(A-T)


Sc = [];
dS = [];
return
if not(matrixFree), % matrix-based
  
  build = isempty(mOld) || isempty(omegaOld) ...
    || isempty(alphaOld) || isempty(lambdaOld) || isempty(muOld) ...
    || length(mOld) ~= length(m) || length(omegaOld) ~= length(omega) ...
    || any(mOld ~= m) || any(omegaOld~=omega) || (alpha ~= alphaOld)...
    || (mu ~= muOld)|| (lambda ~= lambdaOld);
  if build,
    mOld = m; omegaOld = omega; 
    alphaOld = alpha; muOld = mu; lambdaOld = lambda;

    A = getElasticMatrixNodal(omega,m,mu,lambda);
    A = alpha*hd*(A'*A);
  end;
  dS  = uc'*A;
  Sc  = 0.5*dS*uc;
  d2S = A;

else % matrix-free

  d2S.regularizer = regularizer;
  d2S.alpha = alpha;
  d2S.By    = @(u,omega,m) elasticOperatorNodal(u,omega,m,mu,lambda,'By');
  d2S.BTy   = @(u,omega,m) elasticOperatorNodal(u,omega,m,mu,lambda,'BTy');
  d2S.B     = @(omega,m) getElasticMatrixNodal(omega,m,mu,lambda);
  d2S.diag  = @(omega,m) getDiag(omega,m,mu,lambda,alpha);
  
  d2S.d2S   = @(uc,omega,m) ...
    alpha * prod((omega(2:2:end)-omega(1:2:end))./m) * ...
    elasticOperatorNodal(elasticOperatorNodal(uc,omega,m,mu,lambda,'By'),omega,m,mu,lambda,'BTy');
  
  dS        = d2S.d2S(uc,omega,m)';
  Sc        = 0.5*dS*uc;
end

%{
% get diagonal of d2S (interesting in matrix free mode)
function D = getDiag(omega,m,mu,lambda,alpha)
dim   = length(omega)/2;
one   = @(i,j) One(omega,m,i,j);
hd    = prod((omega(2:2:end)-omega(1:2:end))./m);
if dim == 2,
  D   = [reshape( (2*mu+lambda)*one(1,1) + (  mu+lambda)*one(1,2),[],1)
    reshape( (  mu+lambda)*one(2,1) + (2*mu+lambda)*one(2,2),[],1)];
else
  D   = [ ...
    reshape( (2*mu+lambda)*one(1,1)+(mu+lambda)*(one(1,2)+one(1,3)),[],1)
    reshape( (2*mu+lambda)*one(2,2)+(mu+lambda)*(one(2,1)+one(2,3)),[],1)
    reshape( (2*mu+lambda)*one(3,3)+(mu+lambda)*(one(3,1)+one(3,2)),[],1) ];
end;
D = alpha*hd *reshape(D,[],1);

% matrix free implementaion of elastic operator
%
% if strcmp(flag,'By'),
%     compute B  * y
% else
%     compute B' * y
% endif
function By = elasticOperatorNodal(uc,omega,m,mu,lambda,flag)
dim = length(omega)/2; 
h   = (omega(2:2:end)-omega(1:2:end))./m;
hd  = prod((omega(2:2:end)-omega(1:2:end))./m);
a   = sqrt(mu); 
b   = sqrt(mu+lambda);
% nc  = prod(m); % centered
% ns  = prod(ones(length(m),1)*m+eye(length(m)),2); %staggered
nn  = prod(m+1); % nodal 

flag = sprintf('%s-%dD',flag,dim);
switch flag,
  case 'By-2D',
    uc  = reshape(uc,[m+1,2]);
    y1  = squeeze(uc(:,:,1));
    y2  = squeeze(uc(:,:,2));
        
%     p1 = nc;                            %    1:p1   d11 y1
%     p2 = p1+(m(1)+1)*(m(2)-1);          % p1+1:p2   d21 y1
%     p3 = p2+(m(1)-1)*(m(2)+1);          % p2+1:p3   d12 y2
%     p4 = p3+nc;                         % p3+1:p4   d22 y2
%     p5 = p4+nc;                         % p4+1:p5   div y
%     By = zeros(p5,1);

    d1 = @(Y) reshape(Y(2:end,:)-Y(1:end-1,:),[],1)/h(1);
    d2 = @(Y) reshape(Y(:,2:end)-Y(:,1:end-1),[],1)/h(2);
    
    By(1:p4)    = [d1(y1);d2(y1);d1(y2);d2(y2)];
    % div y = d11 y1 + d22 y2
    By(p4+1:p5) = By(1:p1) + By(p3+1:p4);
    By(1:p4)    = a*By(1:p4);
    By(p4+1:p5) = b*By(p4+1:p5);
    
  case 'BTy-2D',
    % the input uc is the result of B*y, therefore it can be splitted into 5
    % parts: d11u1, d21u1,d12u2,d22u2,Divy
    
    p1 = nc;                            %    1:p1   d11 y1
    p2 = p1+(m(1)+1)*(m(2)-1);          % p1+1:p2   d21 y1
    p3 = p2+(m(1)-1)*(m(2)+1);          % p2+1:p3   d12 y2
    p4 = p3+nc;                         % p3+1:p4   d22 y2
    p5 = p4+nc;                         % p4+1:p5   div y
    By = zeros(sum(ns),1);
    
    d1T = @(Y) reshape([-Y(1,:);Y(1:end-1,:)-Y(2:end,:);Y(end,:)],[],1)/h(1);
    d2T = @(Y) reshape([-Y(:,1),Y(:,1:end-1)-Y(:,2:end),Y(:,end)],[],1)/h(2);
    
    By(1:ns(1))   = a*d1T(reshape(uc(   1:p1),m)) ...
      +             a*d2T(reshape(uc(p1+1:p2),m(1)+1,m(2)-1)) ...
      +             b*d1T(reshape(uc(p4+1:p5),m));

    By(1+ns(1):end) = a*d1T(reshape(uc(p2+1:p3),m(1)-1,m(2)+1)) ...
        +             a*d2T(reshape(uc(p3+1:p4),m)) ...
        +             b*d2T(reshape(uc(p4+1:p5),m));   
      
  case 'By-3D'
    y1  = reshape(uc(1:ns(1)),m(1)+1,m(2),m(3));
    y2  = reshape(uc(ns(1)+(1:ns(2))),m(1),m(2)+1,m(3));
    y3  = reshape(uc(ns(1)+ns(2)+(1:ns(3))),m(1),m(2),m(3)+1);
    
    p1  = nc;                            %    1:p1   d11 y1
    p2  = p1+(m(1)+1)*(m(2)-1)*m(3);     % p1+1:p2   d21 y1
    p3  = p2+(m(1)+1)*m(2)*(m(3)-1);     % p2+1:p3   d31 y1
    p4  = p3+(m(1)-1)*(m(2)+1)*m(3);     % p3+1:p4   d12 y2
    p5  = p4+nc;                         % p4+1:p5   d22 y2
    p6  = p5+m(1)*(m(2)+1)*(m(3)-1);     % p5+1:p6   d32 y2
    p7  = p6+(m(1)-1)*m(2)*(m(3)+1);     % p6+1:p7   d13 y3
    p8  = p7+m(1)*(m(2)-1)*(m(3)+1);     % p7+1:p8   d23 y3
    p9  = p8+nc;                         % p8+1:p9   d33 y3
    p10 = p9+nc;                         % p9+1:p10  div y
    By = zeros(p10,1);

    d1 = @(Y) reshape(Y(2:end,:,:)-Y(1:end-1,:,:),[],1)/h(1);
    d2 = @(Y) reshape(Y(:,2:end,:)-Y(:,1:end-1,:),[],1)/h(2);
    d3 = @(Y) reshape(Y(:,:,2:end)-Y(:,:,1:end-1),[],1)/h(3);
    
    By(1:p9)    = [...
      d1(y1);d2(y1);d3(y1);...
      d1(y2);d2(y2);d3(y2);...
      d1(y3);d2(y3);d3(y3) ];
    % div y = d11 y1 + d22 y2 + d33 y3
    By(p9+1:p10) = By(1:p1) + By(p4+1:p5)+ By(p8+1:p9);
    By(   1:p9 ) = a* By(1:p9);
    By(p9+1:p10) = b* By(p9+1:p10);
            
  case 'BTy-3D'
    
    % the input uc is the result of B*y, therefore it can be splitted into 10
    % parts: d11y1, d21y1,d31y1,d12y2,d22y2,d32y2,d13y3,d23y3,u33y3,Divy

    p1  = nc;                            %    1:p1   d11 y1
    p2  = p1+(m(1)+1)*(m(2)-1)*m(3);     % p1+1:p2   d21 y1
    p3  = p2+(m(1)+1)*m(2)*(m(3)-1);     % p2+1:p3   d31 y1
    p4  = p3+(m(1)-1)*(m(2)+1)*m(3);     % p3+1:p4   d12 y2
    p5  = p4+nc;                         % p4+1:p5   d22 y2
    p6  = p5+m(1)*(m(2)+1)*(m(3)-1);     % p5+1:p6   d32 y2
    p7  = p6+(m(1)-1)*m(2)*(m(3)+1);     % p6+1:p7   d13 y3
    p8  = p7+m(1)*(m(2)-1)*(m(3)+1);     % p7+1:p8   d23 y3
    p9  = p8+nc;                         % p8+1:p9   d33 y3
    p10 = p9+nc;                         % p9+1:p10  div y
    By  = zeros(sum(ns),1);
    
    d1T = @(Y) reshape(d1t(Y),[],1)/h(1);
    d2T = @(Y) reshape(d2t(Y),[],1)/h(2);
    d3T = @(Y) reshape(d3t(Y),[],1)/h(3);

     
    By(1:ns(1))   = a*d1T(reshape(uc(   1:p1),m)) ...
      +             a*d2T(reshape(uc(p1+1:p2),m(1)+1,m(2)-1,m(3))) ... 
      +             a*d3T(reshape(uc(p2+1:p3),m(1)+1,m(2),m(3)-1)) ...
      +             b*d1T(reshape(uc(p9+1:p10),m));

    By(ns(1)+(1:ns(2))) = ...
          a*d1T(reshape(uc(p3+1:p4),m(1)-1,m(2)+1,m(3))) ...
      +   a*d2T(reshape(uc(p4+1:p5),m)) ...
      +   a*d3T(reshape(uc(p5+1:p6),m(1),m(2)+1,m(3)-1)) ...
      +   b*d2T(reshape(uc(p9+1:p10),m));    

    By(ns(1)+ns(2)+1:end) = ...
          a*d1T(reshape(uc(p6+1:p7),m(1)-1,m(2),m(3)+1)) ...
      +   a*d2T(reshape(uc(p7+1:p8),m(1),m(2)-1,m(3)+1)) ...
      +   a*d3T(reshape(uc(p8+1:p9),m)) ...
      +   b*d3T(reshape(uc(p9+1:p10),m));    
      
end;

% helper for computation of diag(d2S)
function o = One(omega,m,i,j)
h = (omega(2:2:end)-omega(1:2:end))./m;
m = m + [1:length(m) == i];
o = ones(m)/h(j)^2;
switch j,
  case 1, o(2:end-1,:,:) = 2*o(2:end-1,:,:);
  case 2, o(:,2:end-1,:) = 2*o(:,2:end-1,:);
  case 3, o(:,:,2:end-1) = 2*o(:,:,2:end-1);
end;

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
%}
function runMinimalExample
  % --- 2D --- 
  omega = [0 2 0 3];
  m     = [4,8],[17 19];
  xc    = getNodalGrid(omega,m);
  uc    = xc.^2;%reshape( reshape(xc,[],2)*[0 0;0 1], [], 1);
%   [mbS, mbdS, mbd2S] = feval(mfilename,uc,omega,m,'alpha',1,'mu',1,'lambda',0,'matrixFree',0);
  [mfS, mfdS, mfd2S] = feval(mfilename,uc,omega,m,'alpha',1,'mu',1,'lambda',0,'matrixFree',1);
  return
  
  fprintf('Elastic potential : %f (mb) | %f (mf) \n' , mbS,mfS);
  RE = norm(mbdS-mfdS)/norm(mbdS);
  fprintf('rel.error dS (mb vs. mf) : %e \n' , RE);
  figure(1)
  subplot(1,2,1);
  spy(mbd2S)
  title(sprintf('%s-spy(d2S)-(2D)',mfilename));
  
  
  % --- 3D --- 
  omega = [0 2 0 3 0 2];
  m     = [6 5 7];
  xc    = getStaggeredGrid(omega,m);
  uc    = 1e-2*randn(size(xc));
  [mbS, mbdS, mbd2S] = feval(mfilename,uc,omega,m,'alpha',1,'mu',1,'lambda',0,'matrixFree',0);
  [mfS, mfdS, mfd2S] = feval(mfilename,uc,omega,m,'alpha',1,'mu',1,'lambda',0,'matrixFree',1);
  fprintf('Elastic potential : %f (mb) | %f (mf) \n' , mbS,mfS);
  RE = norm(mbdS-mfdS)/norm(mbdS);
  fprintf('rel.error dS (mb vs. mf) : %e \n' , RE);
  figure(1)
  subplot(1,2,2);
  spy(mbd2S)
  title(sprintf('%s-spy(d2S)-(3D)',mfilename));
  
  

