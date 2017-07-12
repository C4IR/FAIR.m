%==============================================================================
% This code is part of the Matlab-based toolbox LagLDDDM - A Lagrangian Gauss--
% Newton--Krylov Solver for Mass- and Intensity-Preserving Diffeomorphic Image
% Registration
%
% For details and license info see
% - https://github.com/C4IR/LagLDDMM%
%==============================================================================
%
% function [vc,yc,wc,MLhis] = MLLDDMM(MLdata,varargin)
%
% Multi-Level for Large Deformation Diffeomorphic Metric Mapping (LDDMM)
%
% minimizes J^h(y) = D^h(T(y),R) + S^h(y) for h=coarse:fine
%
% Input:
%   ML         coarse to fine representations of the data, see getMultilevel.m
%   varagin    optional parameters, see default parameter
% Output:
%   vc          numerical optimizer (velocity field)
%   yc          optimal transformation
%   wc          optimal parameters for PIR
%   MLhis       iteration history
%   para        additional info from objective function
%
% for level=minLevel:maxLevel,
%   get data(level)
%   if level==minLevel, get wOpt running PIR; end;
%   use pre-registered Yref=trafo(wOpt,X) for regularization
%   if level==minLevel
%     v0 = trafo('get','w0'); % use this as starting guess
%   else
%     get v0 by prologating vopt from coarse to finer level
%   end;
%   get vopt using GaussNewton using v0 as starting guess
% end%For
%
%==============================================================================
%
% optimizer are controlled via [default]
%    PIRobj      [= @PIRobjFctn;         ]
%    PIRpara     [= optPara('PIR-GN');   ]
%    NPIRobj     [= @NPIRobjFctn;         ]
%    NPIRpara    [= optPara('NPIR-GN')   ]
%
% see also Ex_MLLDDMM.m for an example
%==============================================================================

function [vc,yc,wc,MLhis] = MLLDDMM(ML,varargin)

if nargin == 0
  % run minimal example
  help(mfilename);
  Ex_MLLDDMM;
  yc = 'endOfMinimalExample';
  return;
end;

%------------------------------------------------------------------------------
% initialize some default parameter:
% 1.1  initialize output
% 1.2  determine minLevel and maxLevel from the data ML
% 1.3  flag for parametric registration (PIR)
% 2.1  regularization of of PIR objective
% 2.2  regularization of Hessian approximation, H -> H + beta*speye(size(H))
% 2.3  setup stopping value wStop and initial value w0
% 2.4  setup objective in PIR  and optimization parameters
% 3    setup objective in NPIR and optimization parameters
% 4       setup defaults for plots
% 5    replace defaults by user input
% 6    initializes dimension, discretization, grids etc
% 7    convert optimization parameter sets for PIR and NPIR to lists
%------------------------------------------------------------------------------
% 1.1 initialize output
vc          = [];
yc          = [];
wc          = [];
MLhis       = [];

% 1.2  get minLevel and maxLevel from ML
[ML,minLevel,maxLevel] = getMultilevel(ML);

% 1.3  flag for parametric registration (PIR)
parametric  = 1;     % flag: if true, run parametric pre-registration

% 2.1  regularization of PIR objective
%   J_{PIR} = D(T(Y(wc)),Rc) + S(wc) with S_{PIR} = 0.5*(wc-wRef)'*M*(wc-wRef);
wRef        = [];
M           = [];

% 2.2  regularization of Hessian approximation, H -> H + beta*speye(size(H))
beta        = 0;

% 2.3  setup stopping value wStop and initial value w0
wStop       = [];           % global stopping for PIR
w0          = [];           % starting guess for PIR

% 2.4  PIR: default objective function and optimization parameter via optPara('PIR-GN')
PIRobj      = @PIRobjFctn;
PIRpara     = optPara('PIR-GN');

% 3.   NPIR: default objective function and optimization parameter via optPara('NPIR-GN')
NPIRobj     = @LDDMMobjFctn;
NPIRpara    = optPara('NPIR-GN');

% 4       setup defaults for plots
pause       = 0;        % flag for pauses
plots       = 1;        % flag for plots
plotIter    = 0;        % flag for output of iteration history each level
plotMLiter  = 0;        % flag for output of summarized iteration history

omegaV      = [];       % domain for velocity field
N           = 5;        % number of forward euler steps
mV          = @(m) m;   % resolution of velocity field
nt          = regularizer('get','nt');        % number of time discretizations for velocities
vStop       = [];           % global stopping for NPIR
vRef        = [];           % regularization: S(vc-vRef)


% 5    replace defaults by user input
for k=1:2:length(varargin) % overwrites default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

% 6    initializes dimension, discretization, grids etc
dimstr      = @(m) sprintf('[%s]',sprintf(' %d',m));
omega       = ML{end}.omega;  % spatial domain
grid        = regularizer('get','grid');
switch grid
  case 'cell-centered', getGrid = @(m) getCellCenteredGrid(omega,m);
  case 'staggered',     getGrid = @(m) getStaggeredGrid(omega,m);
  case 'nodal',         getGrid = @(m) getNodalGrid(omega,m);
  otherwise, error('nyi');
end;

% 7    convert optimization parameter sets for PIR and NPIR to lists
PO = FAIRcell2struct(PIRpara);
NO = FAIRcell2struct(NPIRpara);
%------------------------------------------------------------------------------

fprintf('%s: MultiLevel Large Deformation Diffeomorphic Metric Mapping\n',mfilename)
fprintf('>> %-20s : %s\n','omega',dimstr(omega));
fprintf('>> %-20s : %s\n','IMAGE MODEL',imgModel);
fprintf('>> %-20s : %s\n','DISTANCE',distance);
fprintf('>> %-20s : %s\n','TRAFO',trafo);
fprintf('>> %-20s : %s, alpha=%s\n','REGULARIZER',...
  regularizer,num2str(regularizer('get','alpha')));
fprintf('>> %-20s : %d\n','PRE-REGISTRATION',parametric);
fprintf('>> %-20s : objective=<%s>, method = <%s>\n',...
  'PIR',func2str(PIRobj),func2str(PIRpara.scheme));
fprintf('>> %-20s : objective=<%s>, method = <%s>\n',...
  'NPIR',func2str(NPIRobj),func2str(NPIRpara.scheme));
fprintf('>> %-20s : minLevel=%d, maxLevel=%d\n','LEVEL',minLevel,maxLevel);
FAIRmessage('#')

tic;
%--------------------------------------------------------------------------
for level=minLevel:maxLevel

  FAIRmessage(sprintf('%s: level %d from %d to %d, %s',...
    mfilename,level,minLevel,maxLevel,dimstr(ML{level}.m)));

  % save(old grid), update(m,grid,data coefficients)
  m     = ML{level}.m;
  xc    = getGrid(m);
  [T,R] = imgModel('coefficients',ML{level}.T,ML{level}.R,omega);
  Rc    = imgModel(R,omega,center(xc,m));
  if strcmp(func2str(NPIRobj),'MPLDDMMobjFctn')
      T = imgModel(T,omega,center(xc,m));
  end

  if level == minLevel    % on coarse level

    if parametric          % perform parametric pre-registration

      % initialize plots for PIR
      FAIRplots('set','mode','PIR','fig',minLevel-1,'plots',plots);
      FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m));

      % ----- call Parametric Image Registration ----------------------------
      xC = getGrid(m);
      PIRfctn = @(wc) PIRobj(T,Rc,omega,m,beta,M,wRef,xC,wc);
      PIRfctn([]); % report status of objective function
      if isempty(w0),    w0    = trafo('w0'); end;
      if isempty(wStop), wStop = w0;          end;

      [wc,hisPIR] = PIRpara.scheme(PIRfctn,w0,'yStop',wStop,PO{:});
      fprintf('wc = \n');
      disp(wc');
      doPause(pause);
      % ----- end PIR --------------------------------------

      if plotIter
        hisPIR.str{1} = 'iteration history for PIR';
        plotIterationHistory(hisPIR,'J',1:4,'fh',100+minLevel-1);
      end;
    else                    % no pre-registration
      wc = [];
    end;

  end;

  % compute xc = trafo(wc,xc)
  if ~isempty(wc) % use yRef(x) = x
    xc = trafo(wc,xc);
  end

  % compute starting guess y0 for this level and stopping value yStop
  if level == minLevel
      if isempty(vRef)
          vRef = getVelocityStartingGuess(omegaV,mV(m),nt);
      end
      v0 = vRef;  % the best known so far
  else
    % prolongate vc (coarse) vRef (current)
    vcOld   = reshape(vc,[],nt+1);
    vRefOld = reshape(vRef,[],nt+1);
    v0 = []; vRef = [];
    for k=1:nt+1
        v0   = [v0 mfPu(vcOld(:,k),omega,mV(m/2))];
        vRef = [vRef mfPu(vRefOld(:,k),omega,mV(m/2))];
    end
    v0 = v0(:); vRef = vRef(:);
 end;
    vStop = 0*getVelocityStartingGuess(omegaV,mV(m),nt);


  % ----- call NPIR -------------------------------------
  % initialize plots for NPIR
  FAIRplots('reset','mode','NPIR','fig',level,'plots',plots);
  FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m));

  NPIRfctn = @(vc) NPIRobj(T,Rc,omega,m,vRef,xc,omegaV,mV(m),N,vc);
  if level == minLevel   % report status of objective function
    NPIRfctn([]);
  end; %

  [vc,his] = NPIRpara.scheme(NPIRfctn,v0,'yStop',vStop,NO{:});
  % ----- end NPIR --------------------------------------

  if plotIter
    his.str{1} = sprintf('iteration history for LDDMM, level=%d',level);
    plotIterationHistory(his,'J',1:4,'fh',100+level);
  end;

  % update iteration history
  if level == minLevel
    MLhis.str = his.str;
    MLhis.his = his.his;
  else
    MLhis.his = [MLhis.his;his.his];
  end;
  doPause(pause);

%--------------------------------------------------------------------------
end;%For level
%--------------------------------------------------------------------------

MLhis.time = toc;

if plotMLiter
  plotMLIterationHistory(MLhis,'fh',[]);
end;

id0 = find(MLhis.his(:,1)==-1,1,'last');
MLhis.reduction = MLhis.his(end,2)/MLhis.his(id0,2);
J = find(MLhis.his(:,1)==-1);
MLhis.iter(minLevel:maxLevel) = MLhis.his([J(2:end)-1;size(MLhis.his,1)],1)';

FAIRmessage([mfilename,' : done !'],'#');

%--------------------------------------------------------------------------

function doPause(p)
if strcmp(p,'on')
  FAIRpause;
elseif p>0
  FAIRpause(p);
end;

%==============================================================================
