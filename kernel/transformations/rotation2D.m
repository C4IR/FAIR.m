%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function [y,dy] = rotation2D(w,x,varargin)
%
% computes y = Q*f(w) and the derivative wrt. w.
% x = reshape(x,[],3); 
% Q = [x(:,1),x(:,2),1, 0      0        0
%      0      0      0, x(:,1),x(:,2,1),1]
% f(w) = codes a rotation around c, dy = Q;
% if no arguments are given, the parameters for the identity map are returned.
%
% see also transformations/contents.m, trafo.m
%==============================================================================

function [y,dy] = rotation2D(w,x,varargin)

% the persitent variable stores the matrix 
% Q(x) = kron( I_2 , [x(:,1),x(:,2),1] );
persistent Q

doDerivative = (nargout>1); % flag for computing the derivative
for k=1:2:length(varargin), % overwrite default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

if nargin==0
    help(mfilename)
    runMinimalExample; 
    return;
else
  y = mfilename('fullfile'); 
  dy = 0;         % parameterization of identity
  if ischar(w),  Q  = []; w = []; end; % reset Q
  if isempty(w), return;          end; 
end;

if isempty(w) || (size(Q,1) ~= numel(x)),
  n = length(x)/2; x = reshape(x,n,2);
  Q = [x,ones(n,1),sparse(n,3);sparse(n,3),x,ones(n,1)];
end;
c = zeros(2,1);
for k=1:2:length(varargin), % overwrites default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

if isempty(Q), feval(mfilename,w,x);  end;

R  = [ cos(w(1)),-sin(w(1));sin(w(1)),cos(w(1))];
g  = (eye(2)-R)*reshape(c,2,1);
f  = [R(1,1);R(1,2);g(1);R(2,1);R(2,2);g(2)];
y  = Q*f;

if ~doDerivative,   return;  end; % no derivative needed     
dR = [-sin(w(1)),-cos(w(1));cos(w(1)),-sin(w(1))];
dg = -dR*reshape(c,[],1);
df = [dR(1,1);dR(1,2);dg(1);dR(2,1);dR(2,2);dg(2)];
dy = Q*df;

%------------------------------------------------------------------------------

function runMinimalExample
fprintf('%s: minimal example\n',mfilename)

omega = [0,10,0,2]; m = [8,9]; 
w = 22/pi;c = (omega(2:2:end)-omega(1:2:end))'/2;
x = getNodalGrid(omega,m);
z = feval(mfilename,w,x,'c',c);

FAIRfigure(1); clf;
plotGrid(x,omega,m,'color','r'); axis image; hold on;
plotGrid(z,omega,m,'color','b'); axis image; hold off;  

fctn = @(w) feval(mfilename,w,x);
w = w + randn(size(w));
checkDerivative(fctn,w,'fig',2);

%==============================================================================
