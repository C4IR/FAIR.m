%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%============================================================================== 
% varargout = viewImage(varargin)
% 
% Main function for image visualization, makes use of persistent parameters
% typical usage:
%   viewPara = {'viewImage','viewImage2D','colormap','gray(256)'};
%   viewImage('reset',viewPara{:});
%   vh = viewImage(T,omega,m);
% see 
%   E2_setupHandsData for a 2D example and 
%   E2_setupBrainData for a 3D example
%
% visualize:
%  
%   visualizes the interpolated data (based on m) on a domain Omega (based on omega)
%
% administration, parameter setting:
%   for resetting, intitializing, setting, updating, clearing, displaying,
%   see dealOptions.m
%
% specific optn:
%   viewImage(T,omega,m,specific{:}), uses specific for this call
%==============================================================================

function varargout = viewImage(varargin)

persistent OPTN 
[method,OPTN,task,stop] = dealOptions(mfilename,OPTN,varargin{:});
if stop,
  varargout{1} = method;
  if nargout > 1, varargout{2} = OPTN;  end;
  return;
end

% update options and split
[method,optn] = dealOptions(mfilename,OPTN,'set',varargin{4:end});
if isempty(method), error('no method specified!'); end;


ih = []; B = [];
% -----------------------------------------------------------------------------
% start to work
Tc     = varargin{1};
omega0 = reshape([0;1]*ones(1,length(size(Tc))),1,[]); 
if nargin<2, omega = omega0;   else omega = varargin{2}; end;
if nargin<3, m     = size(Tc); else m     = varargin{3}; end;

%make optn struct
optn = cell2struct(optn(2:2:end),optn(1:2:end),2);

% shortcut to value of boolean option str
OK = @(str) (isfield(optn,str) && ~isempty(optn.(str)));

% scale image to [0,255]
if OK('scale'), 
  minTc = min(Tc);
  maxTc = max(Tc);
  dG = (maxTc-minTc); dG = dG + 2*(dG == 0);
  Tc = 255/dG*(Tc-minTc);
end;

% invert image
if OK('invert'),
  Tc = max(Tc(:))-Tc; 
  %Tc = 255-max(Tc)+Tc;
end;

% shortcut to value of option str
value = @(str) dealOptions(mfilename,optn,'get',str);

% setup figure
fig = value('fig');
if isempty(fig), 
  fig = gcf;  
elseif fig == 0,
  fig = figure; 
else 
  figure(fig); 
end;

% determine the subfigure
sub = value('sub');
if any(sub ~= 1),
  subplot(sub(1),sub(2),sub(3)); cla;
end;

% set the figure name
figname = value('figname');
if ~isempty(figname),
  set(fig,'numbertitle','off','name',sprintf('[FAIR:%d]: %s',fig,figname));
end;

% splitt options into internal 
internals = {mfilename,'scale','invert','fig','sub','figname'};

% remove internals from list 
for k=1:length(internals),
  if isfield(optn,internals{k}), optn = rmfield(optn,internals{k});  end;
end;

% transfer struct optn back to list
[dummy,optn] = dealOptions(mfilename,optn);

% call the viewer with remaining options
switch nargout,
  case 0, feval(method,Tc,omega,m,optn{:});
  case 1, 
    ih = feval(method,Tc,omega,m,optn{:});
    varargout = {ih}; 
  otherwise, 
    [ih,B]    = feval(method,Tc,omega,m,optn{:});
    varargout = {ih,B}; 
end;
drawnow
%==============================================================================

