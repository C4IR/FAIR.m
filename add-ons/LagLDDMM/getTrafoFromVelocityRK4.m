%==============================================================================
% This code is part of the Matlab-based toolbox 
% LagLDDDM - A Lagrangian Gauss--Newton--Krylov Solver for Mass- and 
%                        Intensity-Preserving Diffeomorphic Image Registration
% 
% For details and license info see 
% - https://github.com/C4IR/FAIR.m/tree/master/add-ons/LagLDDMM
%
% function [yc,dy,para] = getTrafoFromVelocityRK4(vc,yc,varargin)
%
% compute transformation yc by integrating velocity field in time using 
% a fourth-order Runge-Kutta method. Here, vc is assumed to be stationary. 
% For instationary velocities, use getTrafoFromInstatinaryVelocityRK4.m
%
% For more details see Sec. 3 of the paper:
%
% @article{MangRuthotto2017,
%   Title = {A {L}agrangian {G}auss--{N}ewton--{K}rylov solver for mass- and intensity-preserving diffeomorphic image registration},
%   Year = {2017},
%   Journal = {SIAM Journal on Scientific Computing},
%   Author = {A. Mang, L. Ruthotto},
% }
%
% Input:
%
%  vc   - discrete velocity field (nodal, cell-centered, or staggered)
%  yc   - particle positions
%
% Additional REQUIRED Input (provided through varargin)
%
%  omega        - spatial domain (required!)
%  m            - number of cells in each direction (required!)
%
% Optional Input (provided through varargin)
%
%  doDerivative - compute derivative w.r.t. velocity
%  tspan        - time interval (default: [0 1]) 
%  N            - number of time discretization points (nodal) (default: 5)
%  storeInter   - store intermediate transformation (e.g., for visualization)
%
% Output:
%
%  yc           - end point of characteristics
%  dy           - derivative w.r.t. vc
%  para         - info such as CFL
%
% =========================================================================
function [yc,dy,para] = getTrafoFromVelocityRK4(vc,yc,varargin)

if nargin==0
    runMinimalExample
    return;
end
doDerivative = (nargout>1);
omega = [];
m     = [];
tspan = [0 1];
N     = 5;
storeInter = false;
for k=1:2:length(varargin)     % overwrites default parameter
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;
if isempty(omega) || isempty(m)
    error('%s - omega and m must be provided through varargin or trafo(''set'',''omega'',omega,''m'',m')
end
dt   = (tspan(2)-tspan(1))/(N-1); % note dt is negatave when going backwards in time
dy = [];
if doDerivative
    dy = sparse(numel(yc),numel(vc));
end
% compute the CFL number (not actually required for stability in Lagrangian
% methods, but an indicator of how many voxels the particles may move in
% one time step.
h    = (omega(2:2:end)-omega(1:2:end))./m;
% vt   =  sqrt(sum(reshape(vc.^2,[],dim),2));
% CFL  = max(vt)*dt/min(h);
CFL = 0.0;
para = struct('CFL',CFL,'dt',dt,'N',N,'omega',omega,'m',m,'h',h);
if storeInter
    para.YC = zeros(numel(yc),N);
    para.YC(:,1) = yc;
end
for k=1:N-1
    [vi,dvidy] = linearInterGrid(vc,omega,m,yc,'doDerivative',doDerivative);
    if doDerivative
        Ty = getLinearInterGridMatrix(omega,m,yc);
        dyi = Ty + dvidy*dy;
        dytemp = dyi;
    end
    ytemp = vi;
    yi    = yc + .5*dt*vi; 

    [vi,dvidy] = linearInterGrid(vc,omega,m,yi,'doDerivative',doDerivative);
    if doDerivative
        Ty = getLinearInterGridMatrix(omega,m,yi);
        dyi = Ty + dvidy*(dy+.5*dt*dyi);
        dytemp = dytemp + 2*dyi;
    end
    ytemp = ytemp + 2*vi;
    yi    = yc + .5*dt*vi;

    [vi,dvidy] = linearInterGrid(vc,omega,m,yi,'doDerivative',doDerivative);
    if doDerivative
        Ty = getLinearInterGridMatrix(omega,m,yi);
        dyi = Ty + dvidy*(dy+.5*dt*dyi);
        dytemp = dytemp + 2*dyi;
    end
    ytemp = ytemp + 2*vi;
    yi    = yc + dt*vi;
    [vi,dvidy] = linearInterGrid(vc,omega,m,yi,'doDerivative',doDerivative);
    if doDerivative
        Ty = getLinearInterGridMatrix(omega,m,yi);
        dyi = Ty + dvidy*(dy+dt*dyi);
        dytemp = dytemp + dyi;
        dy  = dy + (dt/6)*dytemp;
    end
    ytemp = ytemp + vi;
    yc    = yc + (dt/6)*ytemp;
    if storeInter, para.YC(:,k+1) = yc; end
end




function runMinimalExample

omegaV = [-1 1 -1 1];
omega = .1*omegaV;
m     = [24 24];
tspan     = [2 3];
N     = 2;
yc    = getNodalGrid(omega,m)+.0002;


% test cell-centered grid
regularizer('set','regularizer','mbCurvature','alpha',1)
xc      = reshape(getCellCenteredGrid(omegaV,m),[],2);
vc      = [.3*sign(xc(:,1)).*xc(:,1).^2 sin(pi*xc(:,2))];
fctn    = @(vc) getTrafoFromVelocityRK4(vc,yc,'omega',omegaV,'m',m,'T',tspan,'N',N);
checkDerivative(fctn,vc(:))


% test nodal grid
regularizer('set','regularizer','mbElasticNodal','alpha',1)
xc      = reshape(getNodalGrid(omegaV,m),[],2);
vc      = [.3*sign(xc(:,1)).*xc(:,1).^2 sin(pi*xc(:,2))];
fctn    = @(vc) getTrafoFromVelocityRK4(vc,yc,'omega',omegaV,'m',m,'T',tspan,'N',N);
checkDerivative(fctn,vc(:))

% test staggered grid
regularizer('set','regularizer','mbElastic','alpha',1)
vc  = grid2grid(vc(:),m,'nodal','staggered');
fctn    = @(vc) getTrafoFromVelocityRK4(vc,yc,'omega',omegaV,'m',m,'T',tspan,'N',N);
checkDerivative(fctn,vc(:))

omegaV = [-1 1 -1 1 -1 1];
omega = .1*omegaV;
m     = [16 16 8];
tspan = [4 3];
N     = 10;
yc    = getNodalGrid(omega,m)+.0002;

regularizer('set','regularizer','mbCurvature','alpha',1)
xc      = reshape(getCellCenteredGrid(omegaV,m),[],3);
vc      = [.3*sign(xc(:,1)).*xc(:,1).^2 sin(pi*xc(:,2)) xc(:,3)];
fctn    = @(vc) getTrafoFromVelocityRK4(vc,yc,'omega',omegaV,'m',m,'T',tspan,'N',N);
% [yc,dy] = fctn(vc(:));
checkDerivative(fctn,0*vc(:))

% test nodal grid
regularizer('set','regularizer','mbElasticNodal','alpha',1)
xc      = reshape(getNodalGrid(omegaV,m),[],3);
vc      = [.3*sign(xc(:,1)).*xc(:,1).^2 sin(pi*xc(:,2)) xc(:,3)];
fctn    = @(vc) getTrafoFromVelocityRK4(vc,yc,'omega',omegaV,'m',m,'T',tspan,'N',N);
checkDerivative(fctn,vc(:))

% test staggered grid
regularizer('set','regularizer','mbElastic','alpha',1)
vc  = grid2grid(vc(:),m,'nodal','staggered');
fctn    = @(vc) getTrafoFromVelocityRK4(vc,yc,'omega',omegaV,'m',m,'T',tspan,'N',N);
checkDerivative(fctn,vc(:))

