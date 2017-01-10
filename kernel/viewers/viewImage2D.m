%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% function ih = viewImage2D(T,omega,m,varargin)
%
% THE 2D image viewer, basically calls image.m with the appropriate x1,x2 axis
%
% Input
%   T           discretized image
%   Omega       domain specification
%               Omega = (omega(1),omega(2)) x ... x (omega(2*d-1,omega(2*d)), 
%               spatial dimension d = length(omega/2)
%    m          number of discretization points
%   varargin    optional parameters like {'axis','off'}
%
% Output
%  ih           image handle
% =======================================================================================

function varargout = viewImage2D(T,omega,m,varargin)

if nargin==0
  help(mfilename);
    runMinimalExample; 
    return;
end

if ~exist('m',    'var'), 
    m     = size(T); 
end;
if ~exist('omega','var'), 
    omega = reshape([ones(1,length(m));m],1,[]); 
end;

h  = (omega(2:2:end)-omega(1:2:end))./m;
xi = @(i) (omega(2*i-1)+h(i)/2:h(i):omega(2*i)-h(i)/2)';

ih  = image(xi(1),xi(2),reshape(T,m)'); axis xy image

% the following lines add some nice stuff to the code.
% if varargin = {'title','FAIR','xlabel','x'}
% the code evaluates "title('FAIR');xlabel('x');"
for k=1:2:length(varargin), 
  if ~isempty(varargin{k}), feval(varargin{k},varargin{k+1}); end;
end;

if nargout > 0, varargout = {ih};  end;

%------------------------------------------------------------------------------

function runMinimalExample

setup2DhandData;
FAIRfigure(1); clf;
viewImage2D(dataT,omega,m,'colormap','bone(256)');

%==============================================================================

