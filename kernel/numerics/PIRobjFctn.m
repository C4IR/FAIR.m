%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function  [Jc,para,dJ,H] = PIRobjFctn(T,Rc,omega,m,beta,M,wRef,xc,wc)
%
% Objective Function for Parametric Image Registration
%
% computes J(wc) = D(T(Y(wc)),Rc) + S(wc), where
%
% yc       = y(wc,xc) = trafo(wc,xc)
% Tc       = T(yc)    = imgModel(T,omega,yc), 
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
%
% Output:
%  Jc       current function value J(wc)
%  para     struct {Tc=T(y(wc)), Rc, omega, m, yc=y(wc,xc), Jc}, for plots
%  dJ       gradient of J
%  H        approximation to Hessian of J
%
% see also E6_HNSP_PIR_GN
%==============================================================================

function [Jc,para,dJ,H] = PIRobjFctn(T,Rc,omega,m,beta,M,wRef,xc,wc)

if nargin == 0,

  help(mfilename); 
  runMinimalExample
  Jc = 'endMinimalExample';
  return;

elseif ~exist('wc','var') || isempty(wc),

  % if wc is not an input argument, reports status 
  if nargout == 1, Jc = 'PIR';  return; end;

  % report current settings
  dimstr  = @(m) sprintf('[%s]',sprintf(' %d',m));
  wc      = trafo('w0'); % check for dimensions
  
  FAIRmessage('Parametric Image Registration','-');
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
[yc,dy] = trafo(wc,center(xc,m),'doDerivative',doDerivative);
[Tc,dT] = imgModel(T,omega,yc,'doDerivative',doDerivative);
[Dc,rc,dD,dr,d2psi] = distance(Tc,Rc,omega,m,'doDerivative',doDerivative);

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

% multiply outer and inner derivatives, note: dy might be sparse
if isnumeric(dy)
  
  dD = dD*dT*dy;
  dJ = dD + dS;
  if nargout<4, return; end;

  % approximation to Hessian
  % H = (dy'*(dT'*dT)*dy);
  dr = dr*dT*dy;
  H  = dr'*d2psi*dr + M + beta*speye(length(wc));
  
elseif isstruct(dy),
 
  dD = dy.Qadjoint((dD*dT)')';
  dJ = dD + dS;
  if nargout<4, return; end;
  
  % approximation to Hessian
  H.solver   = 'CG';
  H.operator = @(wc) ...
    dy.Qadjoint(dT'*(dr'*(d2psi* dr*(dT*dy.Q(wc))))) + M*wc + beta*wc;

elseif iscell(dy),

  dD = kronIdy(dD*dT,dy);
  dJ = dD + dS;
  if nargout<4, return; end;

  % approximation to Hessian
  dr = kronIdy(dr*dT,dy);
  H   = dr'*d2psi*dr + M + beta*speye(length(wc));
 
else
  
  error('nyi')

end;

%------------------------------------------------------------------------------

function dd = kronIdy(Df,Dy)
% It is assumed that Q = I_d \otimes Dy{1},
%
% [df1,...dfd]*| Q       | = [df1*Q,...,dfd*Q]
%              |  \      |
%              |        Q|

m   = size(Df,2);
n   = size(Dy{1},1);
switch m/n,
  case 2, dd = [Df(:,1:n)*Dy{:},Df(:,n+1:end)*Dy{:}];
  case 3, dd = [...
      Df(:,1:n)*Dy{:},...
      Df(:,n+1:2*n)*Dy{:},...
      Df(:,2*n+1:end)*Dy{:}];
  otherwise,
    error('nyi')
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
trafo('reset','trafo','rigid2D');
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
OPTpara = FAIRcell2struct(optPara('PIR-GN'));
[wc,his] = GaussNewton(fctn,w0,OPTpara{:}); 
return;

% %------------------------------------------------------------------------------
% %% finally: run the MultiLevel Non-Parametric Image Registration
% [wc,his] = MLPIR(ML);
%==============================================================================
