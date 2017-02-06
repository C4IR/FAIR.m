%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function showResults(ML,yc,varargin)
%
% visualizes registration results
%
% Input:
%   ML          multi-level representation of data
%   yc          transformation
%   varargin    optional parameters like {'isovalue',100}
%==============================================================================

function showResults(ML,yc,varargin)

if nargin==0
    runMinimalExample; return;
end

fig     = [];
figname = sprintf('imgModel=[%s],trafo=[%s],distance=[%s],regularizer=[%s,alpha=%s]',...
  imgModel,trafo,distance,regularizer,num2str(regularizer('get','alpha')));

level  = length(ML);
mode   = 'duplex';
folder = '.';
prefix = mfilename;

for k=1:2:length(varargin), % overwrites defaults
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

omega = ML{level}.omega;
m     = ML{level}.m;
[T,R] = imgModel('coefficients',ML{level}.T,ML{level}.R,omega,'out',0);
yc    = center(yc,m); % make yc cell-centered
xc    = getCellCenteredGrid(omega,m);
R0    = imgModel(R,omega,xc);
T0    = imgModel(T,omega,xc);
T1    = imgModel(T,omega,yc);
d0    = distance(T0,R0,omega,m);
d1    = distance(T1,R0,omega,m);
red   = num2str(100*d1/d0);
D0    = 255 - abs(T0-R0);
D1    = 255 - abs(T1-R0);
FAIRpause(1/100);
FAIRfigure(fig,'figname',figname);

fprintf('reduction %s(yc)/%s(xc)=%s%%\n',distance,distance,red);

if strcmp(mode,'single');
  Name  = @(str) fullfile(folder,sprintf('%s-%s',prefix,str));

  Write = @(I) imwrite(uint8(round(flipud(reshape(I,m)'))),...
    [Name(inputname(1)),'.jpg']);

  viewImage(R0,omega,m); title('R0'); Write(R0); pause(1);
  viewImage(T0,omega,m); title('T0'); Write(T0); pause(1);
  viewImage(D0,omega,m); title('D0'); Write(D0); pause(1);

  clf;   
  viewImage(T0,omega,m); hold on;
  yn = grid2grid(yc,m,'centered','nodal');
  plotGrid(yn,omega,m,'color','w','linewidth',2,...
    'spacing',[max([1,m(1)/32]),max([1,m(2)/32])]);
  axis off;
  FAIRprint(Name('grid'),'obj','gca','folder',[],'pause','off');
  clf;
  viewImage(T1,omega,m); title('T1'); Write(T1); pause(1);
  viewImage(D1,omega,m); title('D1'); Write(D1); pause(1);
  return

else
  subplot(2,3,1); viewImage(R0,omega,m); title('R')
  subplot(2,3,2); viewImage(T0,omega,m); title('T(xc)')
  subplot(2,3,3); viewImage(T1,omega,m); title('T(yc)')
  subplot(2,3,5); viewImage(D0,omega,m);
  title(sprintf('%s(T(xc),R)=%s%%',distance,num2str(100)))
  subplot(2,3,6); viewImage(D1,omega,m);
  title(sprintf('%s(T(yc),R)=%s%%',distance,red))
  subplot(2,3,4); viewImage(T0,omega,m); hold on; title('T(xc) and yc')
  axis(omega)
  plotGrid(yc,omega,m,'spacing',ceil(m/32));
end

%------------------------------------------------------------------------------
function runMinimalExample

close all
% setup data, interpolation, transformation, distance, 
% extract data from a particular level amd compute coefficients for interpolation
setup2DhandData

imgModel('reset','imgModel','splineInter','regularizer','moments','theta',1e-1);
level = 4; 
omega = ML{level}.omega; 
m     = ML{level}.m;
[T,R] = imgModel('coefficients',ML{level}.T,ML{level}.R,omega,'out',0);

distance('reset','distance','SSD');
trafo('reset','trafo','rotation2D','c',(omega(2:2:end)-omega(1:2:end))'/2); 
w0 = trafo('w0'); beta = 0; M =[]; wRef = []; % disable regularization
xc = getCellCenteredGrid(omega,m); 
Rc = imgModel(R,omega,xc);

fctn = @(wc) PIRobjFctn(T,Rc,omega,m,beta,M,wRef,xc,wc);

% optimize
OPTpara = FAIRcell2struct(optPara('PIR-GN'));
wc = GaussNewton(fctn,w0,OPTpara{:},'Plots',0);
yc = trafo(wc,getCellCenteredGrid(ML{end}.omega,ML{end}.m));

showResults(ML,yc);
help(mfilename);
%==============================================================================

