%==============================================================================
% This code is part of the Matlab-based toolbox 
% LagLDDDM - A Lagrangian Gauss--Newton--Krylov Solver for Mass- and 
%                        Intensity-Preserving Diffeomorphic Image Registration
% 
% For details and license info see 
% - https://github.com/C4IR/FAIR.m/tree/master/add-ons/LagLDDMM
%
% function [yc,dy,para] = getTrafoFromInstationaryVelocityRK4(vc,yc,varargin)
%
% compute transformation yc by integrating velocity field in time using a
% Runge-Kutta 4 method with fixed time steps. Velocity is time dependent
% here and assumed to be nodal in time. For stationary velocities, use 
% getTrafoFromVelocityRK4.m
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
%  vc   - discrete velocity field (nodal, cell-centered, or staggered in space / nodal in time)
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
%  N            - number of time steps (default: 5)
%  nt           - number of time points for velocity (default: compute from input)
%  storeInter   - store intermediate transformation (e.g., for visualization)
%
% Output:
%
%  yc           - end point of characteristics
%  dy           - derivative w.r.t. vc
%  para         - struct containing info
%
% =========================================================================
function [yc,dy,para] = getTrafoFromInstationaryVelocityRK4(vc,yc,varargin)

if nargin==0,
    runMinimalExample
    return;
end
doDerivative = (nargout>1);
omega = [];
m     = [];
tspan = [0,1];
N     = 5;
nt    = [];
storeInter = false;
for k=1:2:length(varargin)     % overwrites default parameter
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;
if isempty(omega) || isempty(m)
    error('%s - omega and m must be provided through varargin or trafo(''set'',''omega'',omega,''m'',m')
end
dim = numel(omega)/2;
if isempty(nt) % roughly estimate nt
    nt = round(numel(vc)/(prod(m)*dim))-1;
end
ht    = (tspan(2)-tspan(1))/nt;
tspan = (linspace(tspan(1),tspan(2),N)-min(tspan))/abs(ht)+1; % map to nodal indices
dt    = diff(tspan);
dt    = dt(1);
dy    = [];
if doDerivative
    dy = cell(nt+1,1);
    for j=1:nt+1, dy{j} = sparse(numel(yc),prod(m)*dim); end;
    dyi = dy;
end

para = struct('dt',dt,'N',N,'nt',nt,'omega',omega,'m',m);
if storeInter
    para.YC = zeros(numel(yc),N);
    para.YC(:,1) = yc;
end

vc = reshape(vc,[],nt+1);
for k=1:N-1
    % get v1 = v(yc,t(k));
    tk = tspan(k);
    p  = floor(tk); x = tk-p;  % split x into integer/remainder
    dvi      = zeros(nt+1,1);
    dvi(p)   = 1-x;
    if p<nt+1, dvi(p+1) = x; end
    vt = vc*dvi;
    [vi,dydvi] = linearInterGrid(vt,omega,m,yc,'doDerivative',doDerivative);
    if doDerivative
        Ty = getLinearInterGridMatrix(omega,m,yc);
        for j=1:nt+1
           dyi{j} = dvi(j)*Ty + dydvi*dy{j};
        end
        dytemp  = dyi;
    end
    yi = yc + .5*dt*vi;
	ytemp = vi;
    
    
    % get v2 = v(y1,.5*(t(k+1)+t(k)));
    tk1 = (tspan(k+1)+tspan(k))/2;
    p  = floor(tk1); x = tk1-p;  % split x into integer/remainder
    dvi      = zeros(nt+1,1);
    dvi(p)   = 1-x;
    if p<nt+1, dvi(p+1) = x; end
    vt = vc*dvi;
    [vi,dydvi] = linearInterGrid(vt,omega,m,yi,'doDerivative',doDerivative);
    if doDerivative
        Ty = getLinearInterGridMatrix(omega,m,yi);
        for j=1:nt+1
            dyi{j}    = dvi(j)*Ty + dydvi*(dy{j}+.5*dt*dyi{j});
            dytemp{j} = dytemp{j} + 2*dyi{j};
        end
    end
    yi = yc + .5*dt*vi;
	ytemp = ytemp + 2*vi;
    
    
    % get v3 = v(y2,.5*(t(k+1)+t(k)));
    [vi,dydvi] = linearInterGrid(vt,omega,m,yi,'doDerivative',doDerivative);
    if doDerivative
        Ty = getLinearInterGridMatrix(omega,m,yi);
        for j=1:nt+1
            dyi{j}    = dvi(j)*Ty + dydvi*(dy{j}+.5*dt*dyi{j});
            dytemp{j} = dytemp{j} + 2*dyi{j};
        end
    end
    yi  = yc + dt*vi;
	ytemp = ytemp + 2*vi;

    % get v4 = v(y3,t(k+1)));
    tk1 = tspan(k+1);
    p  = floor(tk1); x = tk1-p;  % split x into integer/remainder
    dvi      = zeros(nt+1,1);
    dvi(p)   = 1-x;
    if p<nt+1, dvi(p+1) = x; end
    vt = vc*dvi;
    [vi,dydvi] = linearInterGrid(vt,omega,m,yi,'doDerivative',doDerivative);
    if doDerivative
        Ty = getLinearInterGridMatrix(omega,m,yi);
        for j=1:nt+1
            dyi{j} = dvi(j)*Ty + dydvi*(dy{j}+dt*dyi{j});
            dytemp{j} = dytemp{j} + dyi{j};
            dy{j}  = dy{j} + (dt/6)*dytemp{j};
        end
    end
    ytemp = ytemp + vi;
    yc    = yc + (dt/6)*ytemp;
    if storeInter, para.YC(:,k+1) = yc; end
end
if doDerivative, dy = horzcat(dy{:}); end;

function runMinimalExample

omegaV = [-1 1 -1 1];
omega = .1*omegaV;
m     = [32 32];
tspan     = [8 0];
N     = 20;
yc    = getNodalGrid(omega,m)+.0002;

% 
% % test cell-centered grid
regularizer('set','regularizer','mbCurvature','alpha',1)
xc      = reshape(getCellCenteredGrid(omegaV,m),[],2);
v0      = [.3*sign(xc(:,1)).*xc(:,1).^2; sin(pi*xc(:,2))];
vc      = v0*[1,1.1,1.5,2];
fctn    = @(vc) getTrafoFromInstationaryVelocityRK4(vc(:),yc,'omega',omegaV,'m',m,'T',tspan,'N',N);
checkDerivative(fctn,vc(:))


