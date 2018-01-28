%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function [y,dy] = splineTransformation2D(w,x,varargin)
%
% computes y = dy*w and returns dy = kron(I2,Q2,Q1), where
% Q{i}(:,1) = spline(x(:,1));
% 
% if no argumanets are given, the parameters for the identity map are returned.
%
% required inputs are: p (number of spline coefficients), omega and m
%
% see also transformations/contents.m, trafo.m 
%==============================================================================


function [y,dy] = splineTransformation2D(w,x,varargin)

% the persitent variable stores the matrix 
% Q(x) = kron( I_2 , Q2, Q1 );
% p is number of spline coefficients, m is size of grid, omega is domain
persistent Q p m omega

% p = []; m = []; omega = []; 
for k=1:2:length(varargin), % overwrite default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

if nargin==0
    help(mfilename)
    runMinimalExample; 
    return;
else
  y = mfilename('fullfile'); 
  dy = zeros(2*prod(p),1);             % parameterization of identity
  if ischar(w),  Q  = []; w = []; end; % reset Q
  if isempty(w), return;          end; 
end;

if isempty(w) || (size(Q,1) ~= numel(x)) || (size(Q,2) ~= numel(w)),
  % it is assumed that x is a cell centered grid, extract xi1 and xi2
  dim = size(omega,2)/2;
  n   = numel(x)/dim;
  if n == prod(m)
    q = m;
  elseif n == prod(m+1)
    q = m+1;
  else
    error('can not handle this grid')
  end;
  
  x  = reshape(x,[q,2]);
  Q1 = getQ1d(omega(1:2),q(1),p(1),x(:,1,1));
  Q2 = getQ1d(omega(3:4),q(2),p(2),x(1,:,2));
  Q  = kron(speye(2),kron(sparse(Q2),sparse(Q1)));
  if nargout == 0, return; end;
end;

y = x(:) + Q*w;
dy = Q;

%------------------------------------------------------------------------------

function Q = getQ1d(omega,m,p,xi)
Q  = zeros(m,p); xi = reshape(xi,[],1);
for j=1:p,
  cj=zeros(p,1); cj(j) = 1;
  Q(:,j) = splineInter(cj,omega,xi);
end;

%------------------------------------------------------------------------------

function runMinimalExample
fprintf('%s: minimal example\n',mfilename)

omega = [0,10,0,8]; m = [8,9]; p = [5,6];
w = zeros([p,2]);  w(3,3,1) = 0.05; w(3,4,2) = -0.1;
x = getCellCenteredGrid(omega,m);
y = feval(mfilename,w(:),x,'omega',omega,'m',m,'p',p,'Q',[]);
FAIRfigure(1); clf;
plotGrid(x,omega,m,'color','r'); axis image; hold on;
plotGrid(y,omega,m,'color','b'); axis image; hold off;  

fctn = @(w) feval(mfilename,w(:),x,'omega',omega,'m',m,'p',p,'Q',[]);
w = w + randn(size(w));
checkDerivative(fctn,w(:),'fig',2);

%==============================================================================
