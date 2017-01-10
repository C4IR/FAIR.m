%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function [LM,fig] = getLandmarks(T,R,omega,m,varargin);
%
% Sets a number q of landmark in 2D template and reference images
%
% Input:
%   T,R         representation of template and reference
%   omega, m    representation of discretization for the domain
%   varargin    optional parameters
%
% Output:
%   LM             q-by-4, the landmarks
%
% see also E5_Hands_TPS for an example
%==============================================================================

function [LM,fig] = getLandmarks(T,R,omega,m,varargin)

if nargin==0
    help(mfilename)
    runMinimalExample; 
    return;
end

LM = [];
for k=1:2:length(varargin),     % overwrites default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

% prepare visualization of T and R
omegaT = omega(1,:);
omegaR = omega(end,:);
xT     = getCellCenteredGrid(omegaT,m);
xR     = getCellCenteredGrid(omegaR,m);
Tc     = imgModel(T,omegaT,xT,'imgModel','linearInter');
Rc     = imgModel(R,omegaR,xR,'imgModel','linearInter');

fig = FAIRfigure; clf; 
if ~isnumeric(fig),
  fig = fig.Number;
end;

set(fig,'numbertitle','off','name',sprintf('[FAIR:%d]',fig));
subplot(1,2,1); viewImage(Tc,omegaT,m); hold on; title(inputname(1));
subplot(1,2,2); viewImage(Rc,omegaR,m); hold on; title(inputname(2));

if isempty(LM) && ~strcmp(FAIRinput,'off'),
    % get a number of landmarks in T
    t = getLM(Tc,omegaT,m,inputname(1),'r',fig,[1,2,1],inf);
    if isempty(t), return;  end;
    
    % get the same number of landmarks in R
    r = getLM(Rc,omegaR,m,inputname(2),'g',fig,[1,2,2],size(t,1));
    if isempty(r) || any(size(t) ~= size(r)), return;  end;
    LM = [t,r];
end;

%------------------------------------------------------------------------------

function LM = getLM(Ic,omega,m,Iname,col,fig,sub,maxNumberLM);

fprintf('set landmarks in %s %s : % s\n',...
  Iname,'[left|middle|right]','[take|abort all|done]');
numberLM = 0;
LM = [];
while 1
  subplot(sub(1),sub(2),sub(3));
  [x,y,b] = ginput(1);
  if isempty(b), b = 3; end;
  switch b,
    case 1, % add new landmark
      numberLM = numberLM + 1;
      LM(numberLM,1:2) = [x,y];
      plot(x,y,[col,'.'],x,y,[col,'o'],'linewidth',2)
      text(x,y,int2str(numberLM),'color',col,'fontsize',20);
    otherwise,
      return;
  end;
  % for the second image
  if numberLM == maxNumberLM, return; end;
end;

%------------------------------------------------------------------------------
function runMinimalExample
setup2DhandData
getLandmarks(dataT,dataR,omega,m);
%==============================================================================
