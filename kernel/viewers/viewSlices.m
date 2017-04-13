%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function varargout = viewSlices(T,omega,m,varargin)
%
% visualizes a 3D image as slice by slice
%
% Input:
%   T           discretized image
%   omega       describing the domain
%    m          number of discretization points
%   varargin    optional parameters like {'axis','off'}
%
% Output:
%   ih          image handle
%   B           the mosaic image
%   frames      number of frames for ij-directions
%==============================================================================

function varargout = viewSlices(T,omega,m,varargin)

if nargin==0
    help(mfilename)
    runMinimalExample; 
    return;
end

% set default parameter
s1 = round(m(1)/2);
s2 = round(m(2)/2);
s3 = round(m(3)/2);
% threshold = min(T(:))+0.1*(max(T(:))-min(T(:)));

for k=1:2:length(varargin), % overwrites default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

T  = reshape(T,m);
h  = (omega(2:2:end)-omega(1:2:end))./m;

% show ortho-slices to first dimension
ih = []; h1 = [];
for i=1:length(s1), 
  x1 = (omega(1) + (s1(i)-.5) * h(1)) * ones(m(2:3));  % fix x1 coordinate
  X  = reshape(getCellCenteredGrid(omega(3:end),m(2:3)),[],2); % 2D grid
  x2 = reshape(X(:,1),m(2:3));
  x3 = reshape(X(:,2),m(2:3));

  Ti = squeeze(T(s1(i),:,:));
  
  sh(i) = surf(x1,x2,x3);
  axis(omega);
  hi = findobj('type','surface');
  h1(end+1) = hi(1);
  set(h1(end),'CData',Ti,'FaceColor','texturemap','edgecolor','none');
  hold on;
end;
ih = [ih;h1]; h1 = [];

% show ortho-slices to second dimension
for i=1:length(s2),
  x2 = (omega(3) + (s2(i)-.5) * h(2)) * ones(m([1 3]));
  X  = reshape(getCellCenteredGrid(omega([1 2 5 6]),m([1 3])),[],2);
  x1 = reshape(X(:,1),m([1 3]));
  x3 = reshape(X(:,2),m([1 3]));
  Ti = squeeze(T(:,s2(i),:));
  
  sh(i) = surf(x1,x2,x3);
  axis(omega);
  hi = findobj('type','surface');
  h1(end+1) = hi(1);
  set(h1(end),'CData',Ti,'FaceColor','texturemap','edgecolor','none');
  hold on;
end;
ih = [ih;h1]; h1 = [];

% show ortho-slices to third dimension
for i=1:length(s3),
  x3 = (omega(5) + (s3(i)-.5) * h(3)) * ones(m(1:2));
  X  = reshape(getCellCenteredGrid(omega(1:4),m(1:2)),[],2);
  x1 = reshape(X(:,1),m(1:2));
  x2 = reshape(X(:,2),m(1:2));
  Ti = squeeze(T(:,:,s3(i)));
  
  sh(i) = surf(x1,x2,x3);
  axis(omega);
  hi = findobj('type','surface');
  h1(end+1) = hi(1);
  set(h1(end),'CData',Ti,'FaceColor','texturemap','edgecolor','none');
  hold on;
end;
hold off;
ih = [ih;h1];

if nargout == 1, varargout = {ih};  end;

%------------------------------------------------------------------------------
function runMinimalExample

load mice3D; 
FAIRfigure(1); clf; hold on
viewSlices(dataT,omega,m);
view(40,45)
%==============================================================================
