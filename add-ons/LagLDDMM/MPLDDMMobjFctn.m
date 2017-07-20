%==============================================================================
% This code is part of the Matlab-based toolbox 
% LagLDDDM - A Lagrangian Gauss--Newton--Krylov Solver for Mass- and 
%                        Intensity-Preserving Diffeomorphic Image Registration
% 
% For details and license info see 
% - https://github.com/C4IR/FAIR.m/tree/master/add-ons/LagLDDMM
%
% function  [Jc,para,dJ,H] = MPLDDMMobjFctn(T,Rc,omega,m,beta,M,wRef,xc,vc)
%
% Objective Function for Mass-preserving LDDMM
%
% computes J(vc) = D(F(yc(vc))*T,Rc) + S(vc - vRef), where
%
% vc       - is a velocity field (stationary and instationary supported)
% F        - is a push forward matrix (paticle-in-cell method)
% yc(vc)   = trafo(wc,xc) is obtained by tracing characteristics
% Tc       = T(yc) = imgModel(T,omega,yc),
% D(Tc,Rc) = distance(Tc,Rc,omega,m)
% S        = regularizer, e.g., S(uc) = 0.5*uc'*B*uc
% vRef     = reference guess for velocities
%
% For more details see the paper:
%
% @article{MangRuthotto2017,
%   Title = {A {L}agrangian {G}auss--{N}ewton--{K}rylov solver for mass- and intensity-preserving diffeomorphic image registration},
%   Year = {2017},
%   Journal = {SIAM Journal on Scientific Computing},
%   Author = {A. Mang, L. Ruthotto},
% }
%
% Input:
%   T       particle mass
%   Rc      sampled reference image
%   omega   spatial domain
%   m       number of discretization points
%   vRef    reference velocity
%   xc      discretization of Omega
%   vc      current velocities
%
% Output:
%  Jc       current function value J(vc)
%  para     struct {Tc=T(y(vc)), Rc, omega, m, yc=y(vc,xc), Jc}, for plots
%  dJ       gradient of J
%  H        approximation to Hessian of J
%
% see also Ex_MPLDDMM.m
%==============================================================================

function [Jc,para,dJ,H] = MPLDDMMobjFctn(T,Rc,omega,m,vRef,xc,omegaV,mV,N,vc)

para = struct([]);
if nargin == 0
    help(mfilename);
    runMinimalExample;
    return;
elseif ~exist('vc','var') || isempty(vc)
    % if wc is not an input argument, reports status
    if nargout == 1, Jc = 'MPLDDMM';  return; end;
    % report current settings
    dimstr  = @(m) sprintf('[%s]',sprintf(' %d',m));
    vc      = trafo('w0');
    fprintf('Mass-preserving LDDMM:\n');
    fprintf('   J(vc) = D(F(yc(vc))*T,Rc) + S(vc-vRef) != min\n');
    fprintf('  %20s : %s\n','m',dimstr(m));
    fprintf('  %20s : %s\n','omega',dimstr(omega));
    fprintf('  %20s : %s\n','IMAGE MODEL',imgModel);
    fprintf('  %20s : %s\n','DISTANCE',distance);
    fprintf('  %20s : %s\n','#timeSteps',num2str(N));
    fprintf('  %20s : %s\n','mV',dimstr(mV));
    fprintf('  %20s : %s\n','omegaV',dimstr(omegaV));
    fprintf('  %20s : %s\n','length(vc)',num2str(length(vc)));
    Jc = vc; % return starting guess
    return;
end;
tspan = [0 1];
dim   = numel(omega)/2;
nt    = round(numel(vc)/(dim*prod(m)))-1;


% do the work ------------------------------------------------------------
matrixFree   = regularizer('get','matrixFree');
doDerivative = (nargout>2);            % flag for necessity of derivatives

% compute transformation, distance, and regularization and combine these
if nt<1
    [yc,dy] = getTrafoFromVelocityRK4(vc,xc,'omega',omegaV,'m',mV,...
                         'N',N,'tspan',tspan,'doDerivative',doDerivative);
else
    [yc,dy] = getTrafoFromInstationaryVelocityRK4(vc,xc,'omega',omegaV,...
                        'm',mV,'N',N,'tspan',tspan,'doDerivative',doDerivative);
end
[I,dI]  = getPICMatrixAnalyticIntegral(omega,m,m,yc,'doDerivative',doDerivative);
Tc      = I*T(:);
if doDerivative, dT      = dI(T(:)); end;

% compute distance
[Dc,~,dD,dres,d2psi] = distance(Tc,Rc,omega,m,'doDerivative',doDerivative);

% compute regularizer
[Sc,dS,d2S] = regularizer(vc-vRef,omegaV,mV,'doDerivative',doDerivative);

Jc = Dc + Sc;

% collect variables for plots
para = struct('Tc',Tc,'Rc',Rc,'omega',omega,'m',m,'yc',yc,'Jc',Jc,'Dc',Dc,'Sc',Sc);

if ~doDerivative, return; end;

dD = dD*dT*dy;
dJ = dD + dS;
if nargout<4, return; end;

dres = dres*dT*dy;
% multiply outer and inner derivatives, note: dy might be sparse
if not(matrixFree)
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
regularizer('reset','regularizer','mfDiffusionST','alpha',[1 1],'nt',2)
v0 = getVelocityStartingGuess(omega,ML{lvl}.m);
xc = getCellCenteredGrid(omega,ML{lvl}.m);
T  = reshape(imgModel(ML{lvl}.T,omega,center(xc,ML{lvl}.m)),ML{lvl}.m);
Rc = imgModel(ML{lvl}.R,omega,center(xc,ML{lvl}.m));

fctn = @(vc) MPLDDMMobjFctn(T,Rc,omega,ML{lvl}.m,v0,xc,omega,ML{lvl}.m,10,vc);
checkDerivative(fctn,v0)


