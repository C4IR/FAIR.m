%==============================================================================
% This code is part of the Matlab-based toolbox LagLDDDM - A Lagrangian Gauss--
% Newton--Krylov Solver for Mass- and Intensity-Preserving Diffeomorphic Image
% Registration
%
% For details and license info see
% - https://github.com/C4IR/LagLDDMM%
%==============================================================================
%
% function  [Jc,para,dJ,H] = LDDMMobjFctn(T,Rc,omega,m,beta,M,wRef,xc,vc)
%
% Objective Function for LDDMM using Lagrangian PDE solver
%
% computes J(vc) = D(T(yc(vc)),Rc) + S(vc - vRef), where
%
% vc       - is a velocity field (stationary or instationary)
% yc(vc)   = trafo(wc,xc) is obtained by tracing characteristics using RK4
%            scheme
% Tc       = T(yc) = imgModel(T,omega,yc),
% D(Tc,Rc) = distance(Tc,Rc,omega,m)
% S        = regularizer, e.g., S(uc) = 0.5*uc'*B*uc
% vRef     = reference guess for velocities
%
% Input:
%   T      - data for template image, Tc = imgModel(T,omega,yc)
%   Rc     - reference image on grid, Rc = imgModel(R,omega,xc)
%   omega  - representation of computational domain
%   m      - discretization size
%   vRef   - reference velocity
%   xc     - discretization of Omega
%   omegaV - computational domain for velocity field (can be larger)
%   mV     - discretization size for velocities
%   N      - number of time steps for RK4 scheme
%   vc     -  current velocities
%
% Output:
%  Jc      - current function value J(vc)
%  para    - struct {Tc=T(y(vc)), Rc, omega, m, yc=y(vc,xc), Jc}, for plots
%  dJ      - gradient of J
%  H       - approximation to Hessian of J
%
% see also ELDDMM_Hand2D
%==============================================================================
function [Jc,para,dJ,H] = LDDMMobjFctn(T,Rc,omega,m,vRef,xc,omegaV,mV,N,vc)

para = struct([]);
if nargin == 0,
  help(mfilename);
  runMinimalExample;
  return;
elseif ~exist('vc','var') || isempty(vc),
  % if wc is not an input argument, reports status
  if nargout == 1, Jc = 'LDDMM';  return; end;
  % report current settings
  dimstr  = @(m) sprintf('[%s]',sprintf(' %d',m));
  vc      = trafo('w0');
  nt      = regularizer('get','nt');
  fprintf('Large Deformation Diffeomorphic Mapping (LDDMM):\n');
  fprintf('   J(vc) = D(T(yc(vc)),Rc) + S(vc-vRef) != min\n');
  fprintf('  %20s : %s\n','m',dimstr(m));
  fprintf('  %20s : %s\n','omega',dimstr(omega));
  fprintf('  %20s : %s\n','IMAGE MODEL',imgModel);
  fprintf('  %20s : %s\n','DISTANCE',distance);
  fprintf('  %20s : %s\n','TRAFO',trafo);
  fprintf('  %20s : %s\n','#timeSteps',num2str(N));
  fprintf('  %20s : %s\n','nt',num2str(nt));
  fprintf('  %20s : %s\n','mV',dimstr(mV));
  fprintf('  %20s : %s\n','omegaV',dimstr(omegaV));
  fprintf('  %20s : %s\n','length(vc)',num2str(length(vc)));
  Jc = vc; % return starting guess
  return;
end;
tspan = [1 0];
dim   = numel(omega)/2;
nt    = round(numel(vc)/(dim*prod(m)))-1;

% do the work ------------------------------------------------------------
matrixFree   = regularizer('get','matrixFree');
doDerivative = (nargout>2);            % flag for necessity of derivatives

% compute transformation, distance, and regularization and combine these
if nt<1
    [yc,dy,pTrafo] = getTrafoFromVelocityRK4(vc,xc,'omega',omegaV,'m',mV,'N',N,'tspan',tspan,'doDerivative',doDerivative);
else
    [yc,dy,pTrafo] = getTrafoFromInstationaryVelocityRK4(vc,xc,'omega',omegaV,...
                        'm',mV,'N',N,'tspan',tspan,'doDerivative',doDerivative);
end
[Tc,dT] = imgModel(T,omega,center(yc,m),'doDerivative',doDerivative);

% compute distance
[Dc,~,dD,dres,d2psi] = distance(Tc,Rc,omega,m,'doDerivative',doDerivative);

% compute regularizer
[Sc,dS,d2S] = regularizer(vc-vRef,omegaV,mV,'doDerivative',doDerivative,'tspan',tspan);

Jc = Dc + Sc;

% collect variables for plots
para = struct('Tc',Tc,'Rc',Rc,'omega',omega,'m',m,'yc',yc,'Jc',Jc,'Sc',Sc,'Dc',Dc,...
    'N',pTrafo.N);

if ~doDerivative, return; end;
dD = dD*dT*dy;
dJ = dD + dS;
if nargout<4, return; end;

% approximation to Hessian
dres = dres*dT*dy;

% multiply outer and inner derivatives, note: dy might be sparse
if not(matrixFree),
  H  = dres'*d2psi*dres + d2S;
else
    % approximation to d2D in matrix free mode
    % d2D   = dr'*d2psi*dr
    % P and P' are operators matrix free
    H.omega     = omegaV;
    H.m         = mV;
    H.nt        = nt;
    H.tspan     = tspan;
    H.d2D.how   = '*dr''*d2psi*dr';
    H.d2D.P     = @(x) x;
    H.d2D.dr    = dres;
    H.d2D.d2psi = d2psi;

    H.d2S = d2S;

end;

function runMinimalExample
setup2DGaussianData;
lvl = 5;
regularizer('reset','regularizer','mbElastic','alpha',1)
v0 = getVelocityStartingGuess(omega,m);
xc = getCellCenteredGrid(omega,MLdata{lvl}.m);
T  = reshape(imgModel(MLdata{lvl}.T,omega,center(xc,MLdata{lvl}.m)),MLdata{lvl}.m);
Rc = imgModel(MLdata{lvl}.R,omega,center(xc,MLdata{lvl}.m));

fctn = @(vc) LDDMMobjFctn(T,Rc,omega,MLdata{lvl}.m,v0,xc,omega,m,10,vc);
checkDerivative(fctn,v0)


