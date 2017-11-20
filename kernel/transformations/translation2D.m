%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% 
% function [y,dy] = translation2D(w,x,varargin)
%
% computes y = Q*w = and the derivative wrt. w.
% x = reshape(x,[],2); 
% Q = [1, 0; 0  1 ];  dy = Q:
% if no arguments are given, the parameters for the identity map are returned.
%
% see also transformations/contents.m, trafo.m 
%==============================================================================

function [y,dy] = translation2D(w,x,varargin)

% the persistent variable Q stores the matrix 
% Q(x) = kron( I_2 , [x(:,1),x(:,2),1] );
persistent Q

if nargin==0
    help(mfilename)
    runMinimalExample; 
    return;
else
  y = mfilename('fullfile'); 
  dy = [0;0];         % parameterization of identity
  if ischar(w),  Q  = []; w = []; end; % reset Q
  if isempty(w), return;          end; 
end;

if isempty(w) || (size(Q,1) ~= numel(x)),
  n = length(x)/2;
  Q = sparse(kron(speye(2),ones(n,1)));
  if nargout == 0, return; end;
end;

y  = x + Q*w;
dy = Q;

%------------------------------------------------------------------------------
function runMinimalExample
fprintf('%s: minimal example\n',mfilename)

omega = [0,10,0,2]; 
m = [8,9]; 
w = (omega(2:2:end)-omega(1:2:end))'/10
x = getNodalGrid(omega,m);
z = feval(mfilename,w,x);

FAIRfigure(1); clf;
plotGrid(x,omega,m,'color','r'); axis image; hold on;
plotGrid(z,omega,m,'color','b'); axis image; hold off;  

fctn = @(w) feval(mfilename,w,x);
w = w + randn(size(w));
checkDerivative(fctn,w,'fig',2);
%==============================================================================
