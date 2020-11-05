%==============================================================================
% This code is part of the Matlab-based toolbox 
% LagLDDDM - A Lagrangian Gauss--Newton--Krylov Solver for Mass- and 
%                        Intensity-Preserving Diffeomorphic Image Registration
% 
% For details and license info see 
% - https://github.com/C4IR/FAIR.m/tree/master/add-ons/LagLDDMM
%
% function [yx,Dyx] = linearInterGrid(yc,omega,m,x)
%
% interpolates grid yc at positions x
%
% Input:
% 
% yc    - grid
% omega - spatial domain
% m     - number of discretization points
% x     - interpolation points
% 
% Output: 
% 
% yx    - y(x), yc evaluated at x
% Dyx   - derivative of y(x) w.r.t. x
%==============================================================================
function [yx, Dyx] = linearInterGrid(yc,omega,m,x,varargin)

if nargin==0
    runMinimalExample;
    return;
end;
inter = @linearInterMex;
doDerivative = (nargout==2);
for k=1:2:length(varargin)    % overwrites default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

dim = length(omega)/2;
m   = m(1:dim);
h   = (omega(2:2:end)-omega(1:2:end))./m;

%% check grid type
cg = numel(yc)/(prod(m)  *dim);
ng = numel(yc)/(prod(m+1)*dim);
sg = numel(yc)/(sum(prod(repmat(m,dim,1)+eye(dim,dim),2)));

if mod(cg,1) == 0
  grid = 'cell-centered';
  yc = reshape(yc,[m,dim]);
  ys = @(d) yc(:,:,d);
  omegaExt = @(p) omega;
  
elseif mod(ng,1) == 0;
  grid = 'nodal';
  yc = reshape(yc,[m+1,dim]);
  ys = @(d) yc(:,:,d);
  omegaExt = @(p) omega + reshape([-h;h]/2,1,[]);
  
elseif mod(sg,1) == 0;
  grid = 'staggered';
  ms  = repmat(m,dim,1)+eye(length(m));
  ns  = prod(ms,2);
  e   = @(p) double((1:length(omega))==p);
  ys  = @(d) reshape(yc(sum(ns(1:d-1))+(1:ns(d))),ms(d,:));
  omegaExt = @(p) omega - e(2*p-1) * h(p)/2 + e(2*p)*h(p)/2;
else 
  error('cannot deal this grid!');
end;

yx  = zeros(numel(x)/dim,dim);
Dyx = [];

for d=1:dim
  [yx(:,d),Dyt] = inter(ys(d),omegaExt(d),x(:),'doDerivative',doDerivative);
  if doDerivative,
    Dyx = [Dyx ; Dyt];
  end
end

yx = yx(:);
          

function runMinimalExample

omega = [0,10,0,8]; m = [8,9]; p = [5,6];
w  = zeros([p,2]);  w(3,3,1) = 0.07; w(3,4,2) = -0.06;
xc = getCellCenteredGrid(omega,m);
xn = getNodalGrid(omega,m);
xs = getStaggeredGrid(omega,m);
yc  = splineTransformation2D(w(:),xn,'omega',omega,'m',m+1,'p',p,'Q',[]);

% some points describing a rectangle
X = [5;6;5;6;2;2;4;4;];

yX = feval(mfilename,yc,omega,m,X);

fctn = @(x) linearInterGrid(yc,omega,m,x);
checkDerivative(fctn,X(:)+randn(size(X(:)))/10);

figure(1); clf;
subplot(1,2,1);
plotGrid(xc,omega,m); hold on; 
plot(X(1:4),X(5:end),'rx');
title('initial grid + points');
subplot(1,2,2);
plotGrid(yc,omega,m); hold on; 
plot(yX(1:4),yX(5:end),'rx');
title('transformed grid + transformed points');