%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% 
% function varargout = FAIRfigure(fig,varargin)
%
% Opens a non-gray (-: FAIR figure at position [position] 
% labels with sets figure number and title [figname]
%
%==============================================================================

function varargout = FAIRfigure(fig,varargin)

caller   = dbstack;        % identify the name of the calling function
caller   = caller(min(length(caller),2)).name;
figname  = caller;
color    = [1,1,0.95];
position = [];

for k=1:2:length(varargin), % overwrite defaults  
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

if nargin == 0,
  fig = [];
end;

if isempty(fig), 
  fig = figureh;  
else
  if ishandle(fig),
    fig = figureh(fig);
  else
    fig = figureh(fig);
  end;
end;

if ~isnumeric(fig),
  fig = fig.Number;
end;

figname  = sprintf('[FAIR:%d] %s',fig,figname);
set(fig,'numbertitle','off','name',figname,'color',color);

if not(isempty(position)) && not(strcmp(position,'default')),
%   keyboard
%     position = FAIRposition('fig',fig,'position',position);
    set(fig,'position',position);
end;
if nargout == 1, varargout = {fig};  end;
%==============================================================================
