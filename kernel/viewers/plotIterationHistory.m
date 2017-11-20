%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% 
% function [fig,ph,th] = plotIterationHistory(his,varargin)
% 
% plots iteration history
% 
% Example: plotIterationHistory(his,'J',[1,2,5],'fig',2);
% Input:
%   his         struct wigth headlines and list
%  varargin     optional parameters, see below
% Output:
%   fig         figure handle
%   ph          plot handle
%   th          text handle
%==============================================================================

function [fig,ph,th] = plotIterationHistory(his,varargin)

if nargin==0
    help(mfilename);
    runMinimalExample; 
    fig =  'endOfMinimalExample';             
    return;
end

fig     = [];
col     = 'krgbcmykrgbcmy';
figname = mfilename;
J       = [];
for k=1:2:length(varargin),
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

if isempty(J), 
    J = 1:length(his.str);  
end;

str = his.str(J); his = his.his(:,J);
dum = abs(max(his,[],1));
dum = dum + (dum==0);

if ~isempty(fig),
  fig = figureh(fig); clf; 
else  
  fig = figureh;      clf; 
end;
if ~isnumeric(fig),
    fig = fig.Number;
end;

set(fig,'numbertitle','off','name',sprintf('[FAIR:%d] %s',fig,figname));
for j=2:size(his,2),
  ph(j-1) = plot(his(2:end,1),his(2:end,j)/dum(j),'color',col(j-1)); hold on;
end;
set(ph,'linewidth',2);
th(1) = title(str{1});
th(2) = xlabel('iter');
th(3) = legend(ph,str{2:end});

%------------------------------------------------------------------------------

function runMinimalExample

% setup data, interpolation, transformation, distance, 
% extract data from a particular level amd compute coefficients for interpolation

setup2DhandData;

level = 4; 
omega = ML{level}.omega; 
m     = ML{level}.m;

imgModel('reset','imgModel','splineInter','regularizer','moments','theta',1e-1);
[T,R] = imgModel('coefficients',ML{level}.T,ML{level}.R,omega,'out',0);

distance('reset','distance','SSD');
trafo('reset','trafo','rotation2D','c',(omega(2:2:end)-omega(1:2:end))'/2); 
w0 = trafo('w0'); beta = 0; M =[]; wRef = []; % disable regularization
xc = getCellCenteredGrid(omega,m); 
Rc = imgModel(R,omega,xc);
fctn = @(wc) PIRobjFctn(T,Rc,omega,m,beta,M,wRef,xc,wc);

% optimize
OPTpara = FAIRcell2struct(optPara('PIR-GN'));    
[~,his] = GaussNewton(fctn,w0,OPTpara{:},'Plots',@(varargin)[]);

plotIterationHistory(his);
%==============================================================================
