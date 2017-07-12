% 2D Example for Parametric Multilevel Mass-Preserving Image Registration using VAMPIRE
% 
% (c) Fabian Gigengack and Lars Ruthotto 2011/02/04, see FAIR.2 and FAIRcopyright.m.
% http://www.uni-muenster.de/EIMI/
% http://www.mic.uni-luebeck.de/
% 
%   - data                 synthetic 2D Gaussian blobs (level 3:8, full
%                          resolution: 256x256)
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             SSD
%   - pre-registration     none
%   - transformation       splineTransformation2Dsparse
%   - regularizer          parametric FAIR regularizer
%   - optimizer            Gauss-Newton with ArmijoBacktrack linesearch

setup2DGaussianData;

% prepare the plot
FAIRplots('clear')
Dshow = @(T,R,omega,m) viewImage2D(128+(T-R)/2,omega,m,'colormap',gray(256));
FAIRplots('set','Dshow',Dshow);

% initialize the regularizer for the non-parametric part
alpha       = 100;
alphaLength = 1;
alphaArea   = 0;
alphaVolume = 1;
[reg,regOptn] = regularizer('reset', 'regularizer', 'mfHyperElastic', ...
    'alpha',alpha, 'alphaLength', alphaLength, 'alphaArea', alphaArea, ...
    'alphaVolume', alphaVolume);

% set the parametric transformation
p = [18 18];
trafo('reset','trafo','splineTransformation2Dsparse','omega',omega,'m',m,'p',p);

% finally: run MLPIR
PIRpara            = optPara('PIR-GN');
PIRpara.solver     = @VAMPIREsolveGN_PCG;

[wc,his] = MLPIR(ML, 'PIRobj', @VAMPIREPIRobjFctn, ...
    'getGrid', @getNodalGrid,'PIRpara',PIRpara);

%% Plot results
yc = trafo(wc,getNodalGrid(omega, m));
% Compute resulting image: dataT(yc) * det(D(yc))
Topt = reshape(linearInter(dataT,omega,center(yc,m)) ...
    .* geometry(yc,m,'Jac','omega',omega), m);
figure;
subplot(2,2,1);
viewImage2D(dataT,omega,m,'colormap', 'gray(256)'); hold on;
plotGrid(center(yc, m), omega, m, 'spacing', [5 5]); axis off; hold off;
title('Template (T) & Grid')
subplot(2,2,2);
viewImage2D(dataR,omega,m,'colormap', 'gray(256)'); axis off;
title('Reference (R)')
subplot(2,2,3);
viewImage2D(Topt,omega,m,'colormap', 'gray(256)'); axis off;
title('VAMPIRE result (Topt)')
subplot(2,2,4);
viewImage2D(abs(Topt-dataR),omega,m,'colormap', 'gray(256)'); axis off;
title('Absolute difference image of R and Topt')