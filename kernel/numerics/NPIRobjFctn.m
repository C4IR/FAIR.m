%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
% ##2
%==============================================================================
% function [Jc,para,dJ,H] = NPIRobjFctn(T,Rc,omega,m,yRef,yc)
%
% This is a template for an objective function to be supplied to an optimization scheme.
% 
% computes J(yc) = D(T(P*yc),Rc) + S(yc-yRef), where
% 
%  Tc       = T(yc) = imgModel(T,omega,P*yc), P a cell-centered grid interpolation of yc
%  D(Tc,Rc) = distance(Tc,Rc,omega,m)
%  S(yc)    = regularizer(yc,omega,m) = 0.5*alpha*hd*norm(B*yc)^2 
%                    
% Modes
%   NPIRobjFctn     - displays the help
%   NPIRobjFctn([]) - reports the current setting of the distance and regularization
%   [Jc,para,dJ,H]  = NPIRobjFctn(T,Rc,omega,m,yRef,yc)
%                   - evaluates the objective function
%                   
% Input
%   T     - data for template image Tc = imgModel(T,omega,yc)
%   Rc    - reference image on grid Rc = imgModel(R,omega,xc)
%   omega - representation of computational domain
%   m     - discretization size
%   yRef  - reference for regularization, Sc = regularizer(yc-yRef,omega,m)
%   yc    - current grid
%   
% Output
%   Jc    - function value
%   para  - parameter for plots
%   dJ    - gradient
%   H     - GaussNewton style approximation to Hessian (either explicit or as operator)
% 
% see also NPIRFGSobjFctn (BFGS version), GaussNewton, E9_Hands_NPIRmb_GN
%==============================================================================

function [Jc,para,dJ,H] = NPIRobjFctn(T,Rc,omega,m,yRef,yc)

persistent P % cell-centered grid interpolation operator

if nargin == 0,

  help(mfilename)
  runMinimalExample;  
  Jc = 'endMinimalExample';
  return
  
elseif ~exist('yc','var') || isempty(yc),

  if nargout == 1, Jc = 'NPIR';  return; end;
  dimstr  = @(m) sprintf('[%s]',sprintf(' %d',m));
  v = @(str) regularizer('get',str); % get regularization configuration like grid

  FAIRmessage('Non-Parametric Image Registration','-'); 
  fprintf('J(yc) = D(T(yc),R) + alpha*S(yc-yReg) != min\n');
  fprintf('%-20s : %s\n','m',dimstr(m));
  fprintf('%-20s : %s\n','omega',dimstr(omega));
  fprintf('%-20s : %s\n','IMAGE MODEL',imgModel);
  fprintf('%-20s : %s\n','DISTANCE',distance);
  fprintf('%-20s : %s\n','REGULARIZER',regularizer);
  fprintf('%-20s : %s\n','alpha,mu,lambda',num2str([v('alpha'),v('mu'),v('lambda')]));
  fprintf('%-20s : %s\n','GRID',v('grid'));
  fprintf('%-20s : %s\n','MATRIX FREE',int2str(v('matrixFree')));
  FAIRmessage('-');

  return;
end;


% do the work ------------------------------------------------------------

% define interpolation for cell-centered grid
P = gridInterpolation(P,omega,m);

doDerivative = (nargout>2);            % flag for necessity of derivatives

% compute interpolated image and derivative, formally: center(yc) = P*yc
[Tc,dT] = imgModel(T,omega,center(yc,m),'doDerivative',doDerivative);

% compute distance measure
[Dc,rc,dD,dres,d2psi] = distance(Tc,Rc,omega,m,'doDerivative',doDerivative);

% compute regularizer
[Sc,dS,d2S] = regularizer(yc-yRef,omega,m,'doDerivative',doDerivative);

% evaluate joint function and return if no derivatives need to be computed
Jc = Dc + Sc;

% store intermediates for outside visualization
para = struct('Tc',Tc,'Rc',Rc,'omega',omega,'m',m,'yc',center(yc,m),'Jc',Jc);

if ~doDerivative, return; end;

if isnumeric(P),
  % matrix based mode
  dr = dres*dT*P;
  dD = dD*dT*P;
  dJ = dD + dS;
  H  = dr'*d2psi*dr + d2S;
else
  % derivatives rather explicit
  dr = dres*dT;
  dD = dD*dT;
  dJ = P(dD')' + dS;
  
  % approximation to d2D in matrix free mode
  % d2D   = P'*dr'*d2psi*dr*P 
  % P and P' are operators matrix free 
  H.omega     = omega;
  H.m         = m;
  H.d2D.how   = 'P''*dr''*d2psi*dr*P';
  H.d2D.P     = P;
  H.d2D.dr    = dr;
  H.d2D.d2psi = d2psi;

  H.d2S = d2S;
end;

%------------------------------------------------------------------------------

function P = gridInterpolation(P,omega,m)

grid = regularizer;
switch grid,
  case 'mbElastic'
    if size(P,1) ~= length(omega)/2*prod(m), % update P
      P = stg2center(m); 
    end;
  case 'mfElastic', P = @(yc) stg2center(yc,m);
  case {'mbCurvature','mbTPS'},  P = 1;         % centered grid - matrix based
  case {'mfCurvature','mfTPS'},  P = @(yc) yc;  % centered grid - matrix free
  case 'mbHyperElastic',
    if size(P,1) ~= length(omega)/2 * prod(m)
      P = nodal2center(m);
    end
  case 'mfHyperElastic',
    P = @(y) nodal2center(y,m);
  otherwise
    grid
    error('can not deak this grid')
end;

%------------------------------------------------------------------------------

function runMinimalExample

setup2DhandData;

% extract data of level 4
level = 4; 
omega = ML{level}.omega; 
m     = ML{level}.m; 

imgModel('reset','imgModel','linearInter'); 
[T,R] = imgModel('coefficients',ML{level}.T,ML{level}.R,omega);
xc    = getCellCenteredGrid(omega,m); 
Rc    = imgModel(R,omega,xc);

% initialize distance measure
distance('reset','distance','SSD');       

% % initialize the transformation and a starting guess
% % trafo('reset','trafo','affine3Dsparse');
% trafo('reset','trafo','rigid2D');
% w0 = trafo('w0');

% initialize the regularization model
regularizer('reset','regularizer','mbElastic','alpha',1000,'mu',1,'lambda',0);

%% build objective function
% note: T  is data for template image
%       Rc is sampled reference image
%       optional Tikhonov-regularization is disabled by setting m = [], wRef = []
%       beta = 0, M = [], wRef = []:  
%       disables additional regularization of Hessian approximation

y0   = getStaggeredGrid(omega,m);
yRef = y0;
fctn = @(yc) NPIRobjFctn(T,Rc,omega,m,yRef,yc); 
fctn([]);   % report status
checkDerivative(fctn,y0+rand(size(y0)));

%% setup plots and initialize
FAIRplots('reset','mode','NPIR-objective','omega',omega,'m',m,'fig',1,'plots',1);
FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m)); 

%% -- solve the optimization problem on one level
OPTpara = FAIRcell2struct(optPara('NPIR-GN'));
[yOpt,his] = GaussNewton(fctn,y0,OPTpara{:},'Plots',@FAIRplots); 

%% finally: run the MultiLevel Non-Parametric Image Registration
% [yc,his] = MLIR(ML);
%==============================================================================
