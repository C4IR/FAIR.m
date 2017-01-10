%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function [y,dy] = rigid3D(w,x,varargin)
%
% computes y = Q*f(w) and the derivative wrt. w.
% x = reshape(x,[],3); 
% Q = [x(:,1),x(:,2),x(:,3),1, 0      0        0      0  0      0        0
%      0      0      0,     0, x(:,1),x(:,2,1),x(:,3),1, 0      0        0  
%      0      0      0,     0, 0      0        0      0, x(:,1),x(:,2,1),x(:,3),1]
% f(w) = combination of Euler angles and translation, dy = Q;
% if no argumanets are given, the parameters for the identity map are returned.
%
% see also transformations/contents.m, trafo.m
%==============================================================================

function [y,dy] = rigid3D(w,x,varargin)

% the persistent variable Q stores the matrix 
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
  dy = zeros(6,1);                     % parameterization of identity
  if ischar(w),  Q  = []; w = []; end; % reset Q
  if isempty(w), return;          end; 
end;

if isempty(w) || (size(Q,1) ~= numel(x)),
  x = reshape(x,[],3);
  Q = sparse(kron(speye(3),[x,ones(size(x,1),1)]));
  if nargout == 0, return; end;
end;

c  = cos(w(1:3)); s  = sin(w(1:3));

R3 = [ c(1),-s(1),    0; s(1), c(1),    0;    0,    0,    1];
R2 = [ c(2),    0, s(2);    0,    1,    0;-s(2),    0, c(2)];
R1 = [    1,    0,    0;    0, c(3),-s(3);    0, s(3), c(3)];
dR3 = [-s(1),-c(1),    0; c(1),-s(1),   0;    0,    0,    0];
dR2 = [-s(2),    0, c(2);    0,    0,   0;-c(2),    0,-s(2)];
dR1 = [    0,    0,    0;    0,-s(3),-c(3);   0, c(3),-s(3)];
f = [R3*R2*R1;reshape(w(4:6),1,3)]; f = f(:);

df = zeros(12,6);
df([1:3,5:7,9:11],:) = [reshape(dR3*R2*R1,9,1),reshape(R3*dR2*R1,9,1),...
  reshape(R3*R2*dR1,9,1),zeros(9,3)];
df([4,8,12],4:6) = eye(3);

y  = Q*f;
if ~doDerivative,   return;  end; % no derivative needed     
dy = Q*df;

%------------------------------------------------------------------------------

function runMinimalExample
help(mfilename);
fprintf('%s: minimal example\n',mfilename)

omega = [0,10,0,8,0,6]; m = [8,7,6];
w = 22/pi;c = (omega(2:2:end)+omega(1:2:end))'/2;
R  = [ cos(w),-sin(w),0;sin(w),cos(w),0;0,0,1]';
g  = (eye(3)-R)*reshape(c,[],1);
w  = [w;0;0;g]
x = getNodalGrid(omega,m);
z = feval(mfilename,w,x);
FAIRfigure(1); clf;
plotGrid(x,omega,m,'color','r'); axis image; hold on;
plotGrid(z,omega,m,'color','b'); axis image; hold off;

fctn = @(w) feval(mfilename,w,x);
w = w + randn(size(w));
checkDerivative(fctn,w,'fig',2);
%==============================================================================

