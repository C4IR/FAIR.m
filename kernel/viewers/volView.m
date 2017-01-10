%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% function varargout = volView(I,omega,m,varargin)
%
% 3D image viewer, basically calls patch.m 
% with appropriate x1,x2,x3 settings
%
% Input:
%   T           discretized image
%   omega       describing the domain
%   m           number of discretization points
%   varargin    optional parameters like {'isovalue',100}
%
% Output:
%   ih          image handle
%==============================================================================

function varargout = volView(I,omega,m,varargin)

if nargin==0
    help(mfilename);
    runMinimalExample; 
    return;
end

% set default parameter
isovalue    = 0;
view        = [-37.5,30];
facecolor   = .75*[1,1,1];
facealpha   = .8;

for k=1:2:length(varargin), % overwrite default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

if nargin < 3, m     = size(I);         end;
if nargin < 2, omega = ones(size(m));   end;

% reshape I
I = reshape(I,m);

% threshold I
I(I<isovalue) = 0;

% prepare geometry
h  = (omega(2:2:end)-omega(1:2:end))./m;
xi = @(i) (omega(2*i-1)+h(i)/2:h(i):omega(2*i)-h(i)/2)';
[x1,x2,x3] = meshgrid(xi(1),xi(2),xi(3));
I = permute(I,[2,1,3]);

cla; 
ih(1) = patch(isosurface(x1,x2,x3,I,isovalue));
ih(2) = patch(isocaps(x1,x2,x3,I,isovalue),'FaceColor','interp');

set(ih(1),...
  'FaceColor',facecolor,...
  'EdgeColor','none',...
  'FaceAlpha',facealpha);

set(ih(2),...
  'EdgeColor','none',...
  'FaceAlpha',facealpha);

isonormals(x1,x2,x3,I,ih(1));
set(ih(1),'FaceColor',[.65 .43 .25],'EdgeColor','none',...
      'AmbientStrength',0.5);
daspect([1 1 1]);

% make view nice
changeView(view)
axis equal vis3d xy
axis(omega);
camlight('right');
%camlight('left');
%camlight(-20,-10);
lighting gouraud
drawnow

if nargout == 1, varargout = {ih};  end;
%------------------------------------------------------------------------------
function changeView(v); view(v);
%------------------------------------------------------------------------------
function runMinimalExample

load mice3D; 
tt = linspace(0,200,21);
FAIRfigure(1); clf;
for i=1:length(tt),
    volView(dataT,omega,m,'isovalue',tt(i),'facecolor',[.85 .4 .25]);
    pause(.1);
end
%==============================================================================
