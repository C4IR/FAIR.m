%==============================================================================
% This code is part of the Finite Element Method app for the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR/FAIRFEM 
%==============================================================================
%
% function [yc,wc,MLhis] = MLIRFEM(ML,varargin)
%
% Multi-Level Image Registration using tetrahedral Finite Elements using
% global refinement
%
% minimizes J^h(y) = D^h(T(y),R) + S^h(y) for h=coarse:fine
% uses PIR [=@GaussNewton] and NPIR [=@GaussNewton]
%               
% Input:
%   ML    coarse to fine representations of the data, see getMultilevel.m
%   varagin     optinonal parameters, see default parameter
% Output:
%   yc          numerical optimizer
%   wc          optimal parameters for PIR
%   MLhis       iteration history
%
% for level=minLevel:maxLevel,
%   get data(level)
%   if level==minLevel, get wOpt running PIR; end;
%   use pre-registered Yref=trafo(wOpt,X) for regularization
%   if level==minLevel
%     y0 = Ystop; % use this as starting guess
%   else
%     get y0 by prologating Yopt from coarse to finer level
%   end;
%   get Yopt running NPIR using y0 as starting guess
% end%For
%==============================================================================

function [yc,wc,MLhis] = MLIRFEM(ML,varargin)

if nargin == 0, help(mfilename); runMinimalExample; return; end;

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
NPIRobj     = @FEMobjFctn;
NPIRpara    = optPara('NPIR-GN');
NPIRpara.Plots = @FAIRplotsFEM;


% 4       setup defaults for plots
pause       = 0;        % flag for pauses
plots       = 1;        % flag for plots
plotIter    = 0;        % flag for output of iteration history each level
plotMLiter  = 0;        % flag for output of summarized iteration history

% 5    replace defaults by user input
for k=1:2:length(varargin), % overwrites default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

% 6    initializes dimension, discretization, grids etc
dimstr      = @(m) sprintf('[%s]',sprintf(' %d',m));
omega       = ML{end}.omega;  % spatial domain
dim         = numel(omega)/2;
xc          = [];                 % current grid

% 7    convert optimization parameter sets for PIR and NPIR to lists
PO = FAIRcell2struct(PIRpara);
NO = FAIRcell2struct(NPIRpara);

getGrid = @(m) getNodalGrid(omega,m);
if dim==2
    getMesh = @TriMesh3;
else
    getMesh = @TetraMesh1;
end
% ---
FAIRmessage(sprintf('%s: MultiLevel Image Registration',mfilename),'#')
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
  xOld  = xc; 
  m     = ML{level}.m; 
  [T,R] = imgModel('coefficients',ML{level}.T,ML{level}.R,omega);
    
  if level == minLevel && parametric, % parametric pre-registration
    % initialize plots
    FAIRplots('reset','mode','PIR','fig',minLevel-1,'plots',plots);
    FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',size(R))); 

    % ----- call Parametric Image Registration ----------------------------
    xC = getGrid(m); 
    Rc    = imgModel(R,omega,center(xC,m));
 
    PIRfctn = @(wc) PIRobj(T,Rc,omega,m,beta,M,wRef,xC,wc); 
    PIRfctn([]); % report status of objective function
    if isempty(w0),    w0    = trafo('w0'); end;
    if isempty(wStop), wStop = w0;          end;
    
    [wc,hisPIR] = PIRpara.scheme(PIRfctn,w0,'yStop',wStop,PO{:});
    
    fprintf('wc = \n'); disp(wc');  doPause(pause);    
    % ----- end PIR --------------------------------------
    if plotIter,                              
      his.str{1} = 'iteration history for PIR'; 
      plotIterationHistory(his,'J',1:4,'fh',100+minLevel-1);
    end;                                      
  elseif level == minLevel, % no pre-registration
     wc = [];                                
  end;                                        
  
  % get triangulation
  Mesh     = getMesh(omega,m);
  xc       = Mesh.xn;
  Rc       = imgModel(R,omega,Mesh.mfPi(xc,'C'));

  % compute yRef = xc or yc = trafo(wc,xc) on the appropriate grid 
  if isempty(wc), % use yRef(x) = x
    yRef = xc;  
  else            % compute yn=y(wc,xNodal) and interpolate on current grid    
    yn   = trafo(wc,xc(:));   
    yRef = reshape(yn,[],dim);  
  end;
  Mesh.xn = reshape(yn,[],dim);
    
  % compute starting guess y0
  if level == minLevel,
    y0 = yRef;  % the best known so far
  else
    % prolongate yc (coarse) y0 (current) 
    y0 = xc(:) + Mesh.mfPuNodal(yc(:) - xOld(:),m/2,'Pu');
    y0 = reshape(y0,[],dim);
  end;

  % ----- call NPIR -------------------------------------
  % S(y) = 0.5*alpha*hd*norm(B*(Y-Yref))^2
  % initialize plots for NPIR
  FAIRplotsFEM('reset','mode','FEM','fig',level,'plots',plots);
  FAIRplotsFEM('init',struct('Tc',T,'Rc',R,'Mesh',Mesh)); 
  NPIRfctn = @(yc) NPIRobj(T,Rc,Mesh,yRef,yc); 
  if level == minLevel, NPIRfctn([]);  end; % report status of objective function
  yStop = yRef;
  [yc,his] = NPIRpara.scheme(NPIRfctn,y0(:),'yStop',yStop(:),NO{:});
  % ----- end NPIR --------------------------------------

  if plotIter,
    his.str{1} = sprintf('iteration history for NPIR, level=%d',level);
    plotIterationHistory(his,'J',1:4,'fh',100+level);
  end;

  % update iteration history
  if level == minLevel,
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
if plotMLiter,
  plotMLIterationHistory(MLhis,'fh',[]);
end;
if isempty(yStop), yStop = xc; end
MLhis.reduction = NPIRfctn(yc)/NPIRfctn(yStop);
J = find(MLhis.his(:,1)==-1); 
MLhis.iter(minLevel:maxLevel) = MLhis.his([J(2:end)-1;size(MLhis.his,1)],1)';

FAIRmessage([mfilename,' : done !']);

function doPause(p)
if strcmp(p,'on'),
  FAIRpause; 
elseif p>0,
  FAIRpause(p);
end;

function runMinimalExample
setup2DhandData

trafo('reset','trafo','affine2D');
distance('reset','distance','SSD');
imgModel('reset','imgModel','splineInter','regularizer','moments','theta',1e-2);
regularizer('reset','regularizer','mbHyperElasticFEM','alpha',1e2);

NPIRpara    = optPara('NPIR-GN');
NPIRpara.Plots = @FAIRplotsFEM;
NPIRpara.lineSearch = @ArmijoDiffeomorphicFEM;
[yc,wc,his] = MLIRFEM(ML,'NPIRpara',NPIRpara,'minLevel',4);

