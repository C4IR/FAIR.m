%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function [Y,Z,c] = YTPS(LM,X,varargin);
%
% Landmark Based Thin-Plate-Spline Transformation
% minimizes sum_{i} { ||Y^i(r_j)-t^i_j|| + theta*S(Y^i) } != min
% where t and r denote the landmarks in template and reference
%
% the solution is explicitly known:
%   y = sum_{j} lambda_j rho(|x-r_j|) + p1(x),
% where rho is the dimension dependent radial basis function and the coefficients 
% a of p1 and the lambda satisfy (A(j,k)=rho(|r_j,r_k|), C = [1,r])
% | A  C | * |lambda| = |t|
% | C' 0 |   |a     |   |0|
%
% Input:
%   LM       position of landmarks [t(:,1),...,t(:,d),r(:,1),...,r(:,d)]
%   X        location of points the TPS has to be evaluated
%   [theta]  optional smoothing parameter, default: theta=0
% Output
%   Y = TPS(X)
%   Z = TPS(r)
%   c = TPS coefficients
%
% see also E5_Hands_TPS for an example
%==============================================================================

function [Y,Z,c] = YTPS(LM,X,varargin)

if nargin==0
    runMinimalExample; return;
end

% initialize parameters
theta = 0;                         
for k=1:2:length(varargin), % overwrite default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

L   = size(LM,1);           % number of landmarks
dim = size(LM,2)/2;         % spacial dimension
t   = LM(:,1:dim);          % LM in template
r   = LM(:,dim+1:2*dim);    % LM in reference

% dimension dependent radial basis function
switch dim,
  case 2,  rho = @(r) r.^2.*log(r+(r==0));
  case 3,  rho = @(r) r;
end;

% setup KKT system and solve it 
A = zeros(L);
for j=1:L,
  for k=1:j-1,
    A(j,k) = rho(norm(r(j,:)-r(k,:)));
  end;
end;
A = A+A';
B = [ones(L,1),r];
KKT = [A+theta*eye(L),B;B',zeros(dim+1)];
c = KKT\[t;zeros(dim+1,dim)];

% evalute Y = y(X)=sum c(k)*rho(|X-r_k|)+p1(X), on grid
% evalute Z = y(r)=sum c(k)*rho(|r-r_k|)+p1(r), on landmarks
% allocate memory
X = reshape(X,[],dim); Y = zeros(size(X)); Z = zeros(L,dim);

% evaluate polynomial part
for i=1:dim,
  Y(:,i) = c(L+1,i)+X*c(L+2:end,i);
  Z(:,i) = c(L+1,i)+r*c(L+2:end,i);
end;

% norm(X-r_j) in compact vectorized version
rad = @(X,j) sqrt(sum((X-ones(size(X))*diag(r(j,:))).^2,2));

% adding radial basis functions one-by-one
for j=1:L,
  for i=1:dim,
    Y(:,i) = Y(:,i) + c(j,i)*rho(rad(X,j));
    Z(:,i) = Z(:,i) + c(j,i)*rho(rad(r,j));
  end;    
end;

% re-organize output
Y = Y(:);
Z = Z(:);
%------------------------------------------------------------------------------
function runMinimalExample
help(mfilename);

setup2DhandData
omegaT = omega(1,:); %#ok<NODEF>
omegaR = omega(end,:);
xT = getCellCenteredGrid(omegaT,m);
xR = getCellCenteredGrid(omegaR,m);
Tc = imgModel(dataT,omegaT,xT);
Rc = imgModel(dataR,omegaR,xR);

theta = 10;
[yc,z] = YTPS(LM(:,1:4),xR,'theta',theta); %#ok<NODEF>
TLM = imgModel(dataT,omegaT,yc);

  
FAIRfigure(1,'figname',mfilename); clf;

subplot(1,3,1); 
  viewImage(Tc,omegaT,m); hold on; axis off;
  ph = plotLM(LM(:,1:2),'numbering','on','color','r');
  set(ph,'linewidth',2,'markersize',20);
  title(sprintf('%s','T&LM'),'fontsize',20);

subplot(1,3,2); 
  viewImage(Rc,omegaR,m); hold on;  axis off;
  ph = plotLM(LM(:,3:4),'numbering','on','color','g','marker','+');
  set(ph,'linewidth',2,'markersize',20);
  title(sprintf('%s','R&LM'),'fontsize',20);

subplot(1,3,3);
  z = reshape(z,[],2);
  viewImage(TLM,omegaR,m); hold on;  axis off;
  ph = plotLM(LM(:,3:4),'numbering','off','color','g','marker','+');
  qh = plotLM(z,'numbering','off','color','m','marker','x');
  rh = plot([LM(:,3) z(:,1)]',[LM(:,4) z(:,2)]','m-','linewidth',3);
  set([ph;qh;rh],'linewidth',2,'markersize',20);
  title(sprintf('T(Y^{TPS},\\theta=%s)&LM',num2str(theta)),'fontsize',30);
%==============================================================================

