%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% function [wOpt,his] = MLPIR(ML,varargin)
%
% Multi-Level Parametric Image Registration
%
% minimizes J^h(w) = D^h(T(Y(w)),R) + S^h(w) for h=coarse:fine
% see PIRobjFctn for a default objective and GaussNewton for a default optimizer
%               
% Input:
%   ML        struct of coarse to fine representations of the data, 
%           see getMultilevel.m
%  varargin optinonal parameters, see default parameter
%
% Output:
%  wOpt     optimal parameter for the parametric part
%  MLhis    iteration history
%
% for level=minLevel:maxLevel,
%   get data(level)
%   if level>minLevel, w0 = wOpt; end;
%   get wOpt running PIR  using w0 as starting guess
% end%For
%==============================================================================

function [wOpt,his] = MLPIR(ML,varargin)

if nargin == 0, 
  help(mfilename); 
  E6_HNSP_MLPIR_SSD_affine2D
  return; 
end;

% setup default parameter for parametric pre-registration
PIRobj      = @PIRobjFctn;
PIRpara     = optPara('PIR-GN');
w0          = trafo('w0');  % starting guess
wStop       = w0;           % starting value (independent of level) used for stopping only

% configure regularizer
wRef        = w0;           % regularization: (w-wRef)'*M*(w-wRef)
M           = [];           %
beta        = 0;            % regularization: H -> H + beta*I
getGrid     = @getCellCenteredGrid;

% setup additional default parameter
pause       = 0;            % flag for pauses
plotIter    = 0;            % flag for output of iteration history each level
plotMLiter  = 0;            % flag for output of summarized iteration history
dimstr      = @(m) sprintf('m = [%s]',sprintf(' %d',m));

[ML,minLevel,maxLevel] = getMultilevel(ML);

for k=1:2:length(varargin),   % overwrites default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

% intialize
omega = ML{end}.omega;
wOpt  = w0;
his   = [];

FAIRmessage(sprintf('%s: MultiLevel Image Registration',mfilename),'#')
fprintf('>> %-20s : %s\n','omega',dimstr(omega));
fprintf('>> %-20s : %s\n','IMAGE MODEL',imgModel);
fprintf('>> %-20s : %s\n','DISTANCE',distance);
fprintf('>> %-20s : %s\n','TRAFO',trafo);
fprintf('>> %-20s :\n','PARARMETRIC-REGISTRATION');
fprintf('>> %-20s : objective=<%s>, method = <%s>\n',...
  'PIR',func2str(PIRobj),func2str(PIRpara.scheme));
fprintf('>> %-20s : minLevel=%d, maxLevel=%d\n','LEVEL',minLevel,maxLevel);
FAIRmessage('#')


FAIRmessage(sprintf('%s, minLevel=%d:maxLevel=%d',...
  'MultiLevel Parametric Image Registration',minLevel,maxLevel));

tic;
% -- for loop over all levels ---------------------------------------------
for level=minLevel:maxLevel;

  FAIRmessage(sprintf('%s: level %d from %d to %d, %s',...
    mfilename,level,minLevel,maxLevel,dimstr(ML{level}.m)));
 
  % get data for current level, compute interpolation coefficients
  m     = ML{level}.m; 
  [T,R] = imgModel('coefficients',ML{level}.T,ML{level}.R,omega);
  
  % update transformation
  trafo('set','omega',omega,'m',m);

  % initialize plots
  PIRpara.Plots('reset','mode','PIR-multi level','fig',level);
  PIRpara.Plots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m)); 
  
  % ----- call PIR ------------------------------------
  xc   = getGrid(omega,m); 
  Rc   = imgModel(R,omega,center(xc,m));
  fctn = @(wc) PIRobj(T,Rc,omega,m,beta,M,wRef,xc,wc);
  if level == minLevel, 
    fctn([]);   % report status
  else
    w0 = wOpt;  % update starting guess
  end; 
  PIRpara.yStop = wStop;
  f = fieldnames(PIRpara);
  v = struct2cell(PIRpara);
  temp = reshape({f{:};v{:}},1,[]);
  [wOpt,hisPIR] = PIRpara.scheme(fctn,w0,temp{:});
  % ----- end PIR --------------------------------------

  if plotIter,                                                
    plotIterationHistory(hisPIR,'J',[1,2,5],'fig',20+level);  
  end;                                                        
  
  % update iteration history
  if level == minLevel,
    his.str = hisPIR.str;
    his.his = hisPIR.his;
  else
    his.his = [his.his;hisPIR.his];
  end;
  doPause(pause)
  
end;%for level
% -- for loop over all levels ---------------------------------------------
his.time = toc;
if plotMLiter,
  plotMLIterationHistory(his,'fig',30);
end;
if isempty(wStop), wStop = w0; end
his.reduction = hisPIR.his(end,2) / (hisPIR.his(1,2)+(hisPIR.his(1,2)==0));
J = find(his.his(:,1)==-1); 
his.iter(minLevel:maxLevel) = his.his([J(2:end)-1;size(his.his,1)],1)';

FAIRmessage([mfilename,' : done !']);

%------------------------------------------------------------------------------

function doPause(p)
if strcmp(p,'on'), 
  pause; 
elseif p>0, 
  pause(p); 
end;
%==============================================================================
