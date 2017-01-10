%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function [fig,pp,tt] = plotMLIterationHistory(his,varargin)
% 
% plots multi-level iteration history
%
% Example: plotMLIterationHistory(his,'fig',3);
%
% Input:
%   his            struct wigth headlines and list
%   varargin    optional parameters, see below
% Output:
%   fig         figure handle
%   ph          plot handle
%   th          text handle
%==============================================================================

function [fig,pp,tt] = plotMLIterationHistory(his,varargin)

if nargin==0
    runMinimalExample; return;
end

J       = [];
fig     = [];
col     = 'krgbcmykrgbcmy';
figname = mfilename;

for k=1:2:length(varargin),
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

if isempty(J), J = 1:length(his.str);  end;
str = his.str(J); his = his.his(:,J);

dum = abs(max(his,[],1));
dum = dum + (dum==0);
if ~isempty(fig),
  fig = figure(fig);  
else  
  fig = figure;     
end;

if ~isnumeric(fig),
  fig = fig.Number;
end;
clf;
if isnumeric(fig),
    fignr = fig;
else
    fignr = fig.Number;
end;
    
set(fig,'numbertitle','off','name',sprintf('[FAIR:%d] %s',fignr,figname));

aa = max(his,[],1);
KK = 1:size(his,1);
JJ = find(his(:,1) < 0);

pp(1)=plot(KK(JJ),his(JJ,2),'s'); hold on;
set(pp(1),'markerfacecolor','b','markersize',10);

pp(2) = plot(KK,his(:,2),'x'); 
set(pp(2),'linewidth',3,'color','b','markersize',10);

pp(3) = plot(KK,his(:,2),'-'); 
set(pp(3),'linewidth',3,'color','b','markersize',10);

% pp=plot(KK(JJ),his(JJ,2),'s',KK,his(:,2),'x'); 

%tt = title('multi-level iteration history: ');
%set(tt,'fontsize',15,'fontweight','bold');
ll = legend(pp(1:2),{'J(y^h_0)','J(y^h_k)'},'location','northeast');
q = max(his(:,2));

for j=1:length(JJ),
  plot(KK(JJ(j))*[1,1],1.2*q*[0,1],'b--'); hold on;
%   tt(1) = text(KK(JJ(j))+0.25,3.5e5,sprintf('%s{-%d}','h=2^',j+2));
  tt(2) = text(KK(JJ(j))+1,q,sprintf('level=%d',j+2));
  set(tt(2),'fontsize',12,'fontweight','normal','rotation',90);
end
set(tt(2),'color','r');
pp(4)=plot(KK(JJ(end):end),his(KK(JJ(end):end),2),'-');
pp(5)=plot(KK(JJ(end)),his(KK(JJ(end)),2),'s');
pp(6)=plot(KK(JJ(end):end),his(KK(JJ(end):end),2),'x');

set(pp(4:6),'linewidth',3,'color','r','markersize',10);
qq = plot(KK(JJ(j))*[1,1],1.2*q*[0,1],'r--');
ll(2) = xlabel('k');
set(ll,'fontsize',15,'fontweight','bold');
set(gca,'fontsize',30,'fontweight','bold');
axis([0,length(his)+1,0,inf])

%------------------------------------------------------------------------------

function runMinimalExample

help(mfilename);

% load data, initialize image viewer, interpolator, transformation, distance
setup2DhandData

imgModel('reset','imgModel','splineInter','regularizer','moments','theta',1e-1);
distance('reset','distance','SSD');
trafo('reset','trafo','affine2D');

% run Multilevel Parametric Image Registration 
[~,his] = MLPIR(ML,'plotIter',0,'plotMLiter',0,'Plots',0);

plotMLIterationHistory(his);
%==============================================================================
