%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================

% function  [Jc,para,dJ,H] = PIRBFGSobjFctn(T,Rc,omega,m,beta,M,wRef,xc,wc)
%
% Objective Function for Parametric Image Registration with lBFGS
%
% computes J(wC) = D(T(Y(wc)),R) + S(wc), where
%
% yc       = y(wc,xc) = trafo(wc,xc)
% Tc       = T(yc) = inter(T,omega,yc), 
% D(Tc,Rc) = distance(Tc,Rc,omega,m)
% S(wc)    = 0.5*(wc-wRef)'*M*(wc-wRef);
%                   
% Input:
%   T       coefficients for template image
%   Rc      sampled reference image
%   omega   spatial domain
%   m       number of discretization points
%   beta    adding beta*I to approximation of Hessian 
%   xc      discretization of Omega
%   wc      current parameters
% Output:
%  Jc       current function value J(wc)
%  para     struct {Tc=T(y(wc)), Rc, omega, m, yc=y(wc,xc), Jc}, for plots
%  dJ       gradient of J
%  H        approximation to Hessian of J
%==============================================================================

function [Jc,para,dJ,H] = PIRBFGSobjFctn(T,Rc,omega,m,beta,M,wRef,xc,wc)

if nargin == 0,

  help(mfilename)
  runMinimalExample;
  Jc = 'endMinimalExample';
  return

elseif ~exist('wc','var') || isempty(wc),

  % if wc is not an input argument, reports status
  if nargout == 1, Jc = 'PIR';  return; end;

  % report current settings
  dimstr  = @(m) sprintf('[%s]',sprintf(' %d',m));
  wc      = trafo('w0');

  FAIRmessage('Parametric Image Registration with BFGS','-');
  fprintf('J(wc)=D(T(y(wc)),R) + (wc-wRef)''*M*(wc-wRef) != min\n');
  fprintf('%-20s : %s\n','m',dimstr(m));
  fprintf('%-20s : %s\n','omega',dimstr(omega));
  fprintf('%-20s : %s\n','IMAGE MODEL',imgModel);
  fprintf('%-20s : %s\n','DISTANCE',distance);
  fprintf('%-20s : %s\n','TRAFO',trafo);
  fprintf('%-20s : %s\n','length(wc)',num2str(length(wc)));
  fprintf('%-20s : %s\n','REGULARIZATION',...
    sprintf('M is %d-by-%d, beta=%s',size(M),num2str(beta))); 
  FAIRmessage('-');

  return;
end;


% do the work ------------------------------------------------------------
doDerivative = (nargout>2);            % flag for necessity of derivatives

% compute transformation, distance, and regularization and combine these
[yc,dy]    = trafo(wc,center(xc,m),'doDerivative',doDerivative);
[Tc,dT]    = imgModel(T,omega,yc,'doDerivative',doDerivative);
[Dc,rc,dD] = distance(Tc,Rc,omega,m,'doDerivative',doDerivative);

% add regularization
if isempty(M) || (all(size(M)==1) && (M==0)),
  Sc = 0;
  dS = 0;
  M  = 0;
else
  dS = (wc-wRef)'*M;
  Sc = 0.5*dS*(wc-wRef);
end;

Jc = Dc + Sc;                           

% collect variables for plots
para = struct('Tc',Tc,'Rc',Rc,'omega',omega,'m',m,'yc',yc,'Jc',Jc);

if ~doDerivative, return; end;
dD = dD * dT;                    % multiply outer and inner derivatives
if size(dD,2) == size(dy,1),     % generic case, dy comes complete
  dJ = dD*dy + dS;     
else                             % tricky case,  dy comes sparse
  n  = size(dy{1},1);
  dJ = reshape(dy{1}'*reshape(dD*dT,n,[]),1,[]) + dS;
%   dr = [dr(:,1:n)*dy{1},dr(:,n+1:2*n)*dy{1},dr(:,2*n+1:3*n)*dy{1}];
end;
if nargout<4, return; end;

% approximation to Hessian
if M == 0,
  H = 1;
else
  H = M;
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

% initialize the transformation and a starting guess
% trafo('reset','trafo','affine3Dsparse');
trafo('reset','trafo','affine2D');
w0 = trafo('w0');

%------------------------------------------------------------------------------
% build objective function
% note: T  is data for template image
%       Rc is sampled reference image
%       optional Tikhonov-regularization is disabled by setting m = [], wRef = []
%       beta = 0, M = [], wRef = []:  
%       disables additional regularization of Hessian approximation
beta = 0; M = []; wRef = [];
fctn = @(wc) PIRobjFctn(T,Rc,omega,m,beta,M,wRef,xc,wc); 
fctn([]);   % report status
checkDerivative(fctn,w0+rand(size(w0)));

% setup plots and initialize
FAIRplots('reset','mode','PIR-objective','omega',omega,'m',m,'fig',1,'plots',1);
FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m)); 

%% -- solve the optimization problem on one level
OPTpara = FAIRcell2struct(optPara('lBFGS','solver',''));
[wc,his] = lBFGS(fctn,w0,OPTpara{:}); 
return;

% %------------------------------------------------------------------------------
% %% finally: run the MultiLevel Non-Parametric Image Registration
% [wc,his] = MLPIR(ML);
%==============================================================================
