%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function [y,dy] = rigid2Dsparse(w,x,varargin)
%
% computes y = Q*f(w) and the derivative wrt. w.
% Q = kron( I2, [x(:,1),x(:,2),1] ); dy = Q;
% f(w) = [cos(w(1));-sin(w(1));w(2);sin(w(1));cos(w(1));w(3)]
% if no arguments are given, the parameters for the identity map are returned.
%
% see also transformations/contents.m, trafo.m
%==============================================================================


function [y,dy] = rigid2Dsparse(w,x,varargin)

% the persitent variable stores the matrix 
% Q(x) = kron( I_2 , [x(:,1),x(:,2),1] );
persistent Q

doDerivative = (nargout>1); % flag for computing the derivative
for k=1:2:length(varargin), % overwrite default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

if nargin == 0, 
  runMinimalExample;
  return;
else
  y = mfilename('fullfile'); 
  dy = [0;0;0];         % parameterization of identity
  if ischar(w),  Q  = []; w = []; end; % reset Q
  if isempty(w), return;          end; 
end;

% test for need of updating Q
OK = ~any(isempty(Q));
if OK,
  m1 = size(Q{1});  
  OK = (2*m1(1) == size(x,1)) && (2*m1(2) == numel(w));
end;

if ~OK,
  n = length(x)/2; 
  Q = {[reshape(x,[],2),ones(n,1)]};
  if nargout == 0, return; end;
end;

f = [cos(w(1));-sin(w(1));w(2);sin(w(1));cos(w(1));w(3)];
f = reshape(f,3,2);
y = reshape(Q{1}*f,[],1);

if ~doDerivative,   return;  end; % no derivative needed
df = [-sin(w(1)),0,0;-cos(w(1)),0,0;0,1,0;...
       cos(w(1)),0,0;-sin(w(1)),0,0;0,0,1];
df = reshape(df,3,6);
dy = reshape(Q{1}*df,[],3);
%------------------------------------------------------------------------------

function runMinimalExample
help(mfilename);
fprintf('%s: minimal example\n',mfilename)

omega = [5,10,1,2]; m = [8,9]; 
w = 22/pi;c = (omega(2:2:end)+omega(1:2:end))'/2;
R = [ cos(w),-sin(w);sin(w),cos(w)];
g = (eye(2)-R)*reshape(c,2,1);
w = [w;g]
x = getNodalGrid(omega,m);
z = feval(mfilename,w,x);
t = rigid2D(w,x);
n = norm(z-t);
fprintf('diff between rigid2D and rigid2Dsparse is %s\n',num2str(n));

assert(n<1e-13,'difference to non-sparse version too big')

FAIRfigure(1); clf;
plotGrid(x,omega,m,'color','r'); axis image; hold on;
plotGrid(z,omega,m,'color','b'); axis image; hold off;  
%==============================================================================

