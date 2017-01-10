%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% 
% function varargout = viewImage2Dsc(T,omega,m,varargin)
%
% 2D image viewer, basically calls imagesc.m with the approprixate x1,x2
%
% Input:
%   T           discretized image
%   omega       describing the domain
%    m          number of discretization points
%   varargin    optional parameters like {'axis','off'}
%
% Output:
%  ih           image handle
%==============================================================================

function varargout = viewImage2Dsc(T,omega,m,varargin)

if nargin==0
    help(mfilename)
    runMinimalExample; 
    return;
end

h  = (omega(2:2:end)-omega(1:2:end))./m;
xi = @(i) (omega(2*i-1)+h(i)/2:h(i):omega(2*i)-h(i)/2)';
ih = imagesc(xi(1),xi(2),reshape(T,m)'); axis xy image

% the following lines add some nice stuff to the code.
% if varargin = {'title','FAIR','xlabel','x'}
% the code evaluates "title('FAIR');xlabel('x');"
for k=1:2:length(varargin), 
  if ~isempty(varargin{k}), feval(varargin{k},varargin{k+1}); end;
end;
if nargout == 1, varargout = {ih};  end;

%------------------------------------------------------------------------------

function runMinimalExample
setup2DhandData; 
FAIRfigure(2); clf;
viewImage2Dsc(dataT,omega,m,'colormap','bone(256)');

%==============================================================================
