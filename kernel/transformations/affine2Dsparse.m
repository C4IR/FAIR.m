%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function [y,dy] = affine2Dsparse(w,x,varargin)
%
% computes y = kron(I2,Q)*w = [w(1),w(2);w(4),w(5)] * x + w([3,6]) 
% and the derivative wrt. w.
% Q = {[x(:,1),x(:,2),1]}, dy = Q;
% if no argumanets are given, the parameters for the identity map are returned.
%
% see also transformations/contents.m, trafo.m
%==============================================================================


function [y,dy] = affine2Dsparse(w,x,varargin)

% the persitent variable stores the matrix 
% Q(x) = kron( I_2 , [x(:,1),x(:,2),1] );
persistent Q

if nargin == 0, 
  runMinimalExample;
  return;
else
  y = mfilename('fullfile'); 
  dy = [1;0;0;0;1;0];         % parameterization of identity
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

w  = reshape(w,3,2);
% mimicing Qfull*w as [Q*w(:,1),Q*w(:,2)]
y  = [Q{1}*w(:,1);Q{1}*w(:,2)];
dy = Q;

%------------------------------------------------------------------------------
function runMinimalExample
help(mfilename);
fprintf('%s: minimal example\n',mfilename)

omega = [0,10,0,2]; m = [8,9];
w = 22/pi;c = (omega(2:2:end)-omega(1:2:end))'/2;
R = [ cos(w),-sin(w);sin(w),cos(w)];
g = (eye(2)-R)*reshape(c,2,1);
w = [R(1,1);R(1,2);g(1);R(2,1);R(2,2);g(2)]
x = getNodalGrid(omega,m);
z = feval(mfilename,w,x);
t = affine2D(w,x);
n = norm(z-t);
fprintf('diff between affine2D and affine2Dsparse is %s\n',num2str(n));

assert(n<1e-13,'difference to non-sparse version too big')
FAIRfigure(1); clf;
plotGrid(x,omega,m,'color','r'); axis image; hold on;
plotGrid(z,omega,m,'color','b'); axis image; hold off;
%==============================================================================


