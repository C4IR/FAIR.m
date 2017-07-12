% =========================================================================
% ##1
%
% 2D Multilevel Mass-Preserving LDDMM Example using nonstationary velocity
% field and diffusion regularizer.
%
% =========================================================================
close all; clc; clear all;
setup2DGaussianData;
alpha = [1e3 0];
lvl = 5;
pad = 0.4;
mV = @(m) m;
nt = 1;
N  = 2;
minLevel = 5;
maxLevel = 8;
imgModel('reset','imgModel','linearInterMex')
% 1) setup grid for velocities (padded)
omegaV = omega; omegaV(1:2:end) = omegaV(1:2:end)-pad;  omegaV(2:2:end) = omega(2:2:end)+pad;

% 2) setup regularizer (and decide for stationary or nonstationary velocity)
regularizer('reset','regularizer','mfDiffusionST','alpha',alpha,'nt',nt,'HessianShift',1e-2); % nonstationary velocity

%%
plots = 0;
NPIRpara    = optPara('NPIR-GN');
NPIRpara.maxIter = 40;
NPIRpara.scheme = @GaussNewtonLDDMM;

[vc,~,wc,his] = MLLDDMM(ML,'minLevel',minLevel,'maxLevel',maxLevel,'omegaV',omegaV,...
    'mV',mV,'N',N,'parametric',0,'plots',plots,'NPIRpara',NPIRpara,'NPIRobj',@MPLDDMMobjFctn);
        yInv = getTrafoFromInstationaryVelocityRK4(vc,getNodalGrid(omega,m),'omega',omegaV,'m',m,'N',N,'nt',nt,'tspan',[1 0]);

%%
yc = getTrafoFromInstationaryVelocityRK4(vc,getNodalGrid(omega,m),'omega',omegaV,'m',m,'nt',nt,'tspan',[1,0],'N',N);
Jac = geometry(yInv,m,'Jac','omega',omega);
Topt = linearInterMex(dataT,omega,center(yInv,m)) .* Jac;
D0  = distance(dataT(:),dataR(:),omega,m);
DOpt = distance(Topt(:),dataR(:),omega,m);

fig = figure(); clf;
fig.Name = sprintf('Results for %s',mfilename);

subplot(2,3,1);
viewImage(dataR,omega,m);
title('reference');

subplot(2,3,4);
viewImage(dataT,omega,m);
hold on;
plotGrid(yc,omega,m,'spacing',4)
title('template');

subplot(2,3,2);
viewImage(Topt,omega,m);
title('T(yc)')

subplot(2,3,3);
viewImage(dataT(:)-dataR(:),omega,m);
title('init. residual, SSD=100%');

subplot(2,3,5);
viewImage2Dsc(Jac,omega,m);
title(sprintf('Jac, min=%1.2f max=%1.2f',min(Jac(:)),max(Jac(:))));

subplot(2,3,6);
viewImage(Topt(:)-dataR(:),omega,m);
title(sprintf('opt residual, SSD=%1.2f%%',100*DOpt/D0));