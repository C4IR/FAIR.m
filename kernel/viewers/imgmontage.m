%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% function [ih,B,frames] = imgmontage(T,omega,m,varargin)
%
% visualizes a 3D image as montage
%
% Input:
%   I           discretized image
%   omega       domain specification
%    m          number of discretization points
%   varargin    optional parameters like {'axis','off'}
%
% Output:
%   ih          image handle
%   B           the mosaiced image
%   frames      number of frames for ij-directions
%==============================================================================

function [ih,B,frames] = imgmontage(T,omega,m,varargin)

if nargin==0
    help(mfilename)
    runMinimalExample; 
    ih = 'endOfMinimalExample';
    return;
end

% setup default parmeters
threshold = 0;
framesx   = [];
framesy   = [];
viewer2D  = 'imagesc';
direction = 'xyz';
numbering = 'off';
colormap  = [];
title     = [];
xlabel    = [];
ylabel    = [];
color     = 'b';
linewidth = 2;
fontsize  = 15;
textcolor = 'b';
plots     = 1;
coloraxis = [];

if ~exist('m','var'),
  m = size(T);
end;
if ~exist('omega','var'),
  omega = reshape([ones(1,length(m));m],1,[]);
end;

for k=1:2:length(varargin), % overwrites default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;
ih = []; B = []; frames = [];

% threshold image
if threshold ~= 0, T(T<threshold) = 0; end;

T = flipdim(permute(reshape(T,m),[2,1,3]),1);     % ! note: axis image is to be used!

% reorganize the data to any desired directional view 
switch direction,
  case  'xyz', 
  case '-xyz', T = flipdim(T,1);
  case '+xyz', T = flipdim(T,2);
  case  'xzy', T = permute(T,[1,3,2]);
  case '-xzy', T = flipdim(permute(T,[1,3,2]),1);
  case  'yxz', T = permute(T,[2,1,3]);
  case '-yxz', T = flipdim(permute(T,[2,1,3]),1);
  case '+yxz', T = flipdim(permute(T,[2,1,3]),2);
  case  'yzx', T = permute(T,[2,3,1]);
  case '-yzx', T = flipdim(permute(T,[2,3,1]),1);
  case  'zxy', T = permute(T,[3,1,2]);
  case '-zxy', T = flipdim(permute(T,[3,1,2]),1);
  case  'zyx', T = permute(T,[3,2,1]);
  case '-zyx', T = flipdim(permute(T,[3,2,1]),1);
  otherwise, error(['direction=',direction]);
end;

m = size(T);

% compute the frames in x-y directions
if isempty(framesx) & isempty(framesy),
  framesy = ceil(sqrt(m(3))); 
  framesx = ceil(m(3)/framesy);
elseif isempty(framesx),  
  framesx = ceil(m(3)/framesy);
else                      
  framesy = ceil(m(3)/framesx);
end;
% fprintf('framesx=%d,framesy=%d\n',framesx,framesy);
frames = [framesx,framesy];

% allocate memory for mosaic image B and fill it
B = zeros(m(1)*framesx,m(2)*framesy);
px = 0; py = 0;
for i=1:m(3),
  B(px+(1:m(1)),py+(1:m(2))) = squeeze(T(:,:,i));
  py = py + m(2);
  if ~rem(i,framesy), px = px + m(1); py = 0; end;
end;

% return the mosaiced image and do no plots
if plots==0, 
    return; 
end;                

% initialize the figure and do the plots
fig = gcf; cla;
ih  = feval(viewer2D,B);
if not(isempty(coloraxis));    caxis(coloraxis);    end;
axis ij image off
axis(0.5+[0 size(B,2),0,size(B,1)])

% plot lines
hold on
ph = plot(0.5+[0;size(B,2)]*ones(1,framesx+1),0.5+[1;1]*(0:m(1):size(B,1)),'-',...
     0.5+[1;1]*(0:m(2):size(B,2)),0.5+[0;size(B,1)]*ones(1,framesy+1),'-');
set(ph,'color',color,'linewidth',linewidth);
hold off;

% add some stuff to the figure
annotate('colormap',colormap,'title',title,...
  'xlabel',xlabel,'ylabel',ylabel,'zlabel',ylabel);

% add numbers to the slices
if strcmp(numbering,'on'),
  px = 0; py = 0;
  for i=1:m(3),
    text(py+0.1*m(2),px+0.8*m(1),['#',int2str(i)],...
      'fontsize',fontsize,'color',textcolor);
    py = py + m(2);
    if ~rem(i,framesy), px = px + m(1); py = 0; end;
  end;  
end;

%------------------------------------------------------------------------------

function annotate(varargin)
% the following lines add some nice stuff to the code.
% if varargin = {'title','FAIR','xlabel','x'}
% the code evaluates "title('FAIR');xlabel('x');"
for k=1:2:length(varargin), 
  if ~isempty(varargin{k+1}), feval(varargin{k},varargin{k+1}); end;
end;

if nargout == 1, 
  varargout = {ih};
elseif nargout == 2,
  varargout = {ih,B};  
elseif nargout == 3,
  varargout = {ih,B,frames};  
end

%------------------------------------------------------------------------------
function runMinimalExample
load mice3D; 
FAIRfigure(3); clf;
imgmontage(dataT,omega,m,'numbering','on','textcolor','w');
%==============================================================================
