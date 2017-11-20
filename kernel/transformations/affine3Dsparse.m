%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function [y,dy] = affine3Dsparse(w,x,varargin)
%
% computes y = kron(I3,Q)*w and the derivative wrt. w.
% x = reshape(x,[],3); 
% Q = [x(:,1),x(:,2),x(:,3),1], dy = {Q};
% if no arguments are given, the parameters for the identity map are returned.
%
% see also transformations/contents.m, trafo.m
%==============================================================================

function [y,dy] = affine3Dsparse(w,x,varargin)

% the persitent variable stores the matrix 
% Q(x) = kron( I_2 , [x(:,1),x(:,2),1] );
persistent Q

if nargin==0
    help(mfilename)
    runMinimalExample; 
    return;
else
  y = mfilename('fullfile'); 
  dy = reshape(eye(4,3),[],1);         % parameterization of identity
  if ischar(w),  Q  = []; w = []; end; % reset Q
  if isempty(w), return;          end; 
end;

% test for need of updating Q
OK = ~any(isempty(Q));
if OK, 
  m1 = size(Q{1});  
  OK = (3*m1(1) == size(x,1)) && (3*m1(2) == numel(w));
end;

if ~OK,
  n = length(x)/3; 
  Q = {[reshape(x,[],3),ones(n,1)]};
  if nargout == 0, return; end;
end;

w  = reshape(w,4,3);
% mimicing Qfull*w as [Q*w(:,1),Q*w(:,2),Q*w(:,3)]
y  = [Q{1}*w(:,1);Q{1}*w(:,2);Q{1}*w(:,3)];
dy = Q;

%------------------------------------------------------------------------------

function runMinimalExample
fprintf('%s: minimal example\n',mfilename)

omega = [0,10,0,8,0,6]; m = [8,7,6];
w = 22/pi;c = (omega(2:2:end)-omega(1:2:end))'/2;
R = [ cos(w),-sin(w),0;sin(w),cos(w),0;0,0,1];
g = (eye(3)-R)*reshape(c,[],1);
w = reshape([R,g]',[],1)
x = getNodalGrid(omega,m);
z = feval(mfilename,w,x);
t = affine3D(w,x);
n = norm(z-t);
fprintf('diff between affine3D and affine3Dsparse is %s\n',num2str(n));

assert(n<1e-13,'difference to non-sparse version too big')

FAIRfigure(1); clf;
plotGrid(x,omega,m,'color','r'); axis image; hold on;
plotGrid(z,omega,m,'color','b'); axis image; hold off;   
%==============================================================================
