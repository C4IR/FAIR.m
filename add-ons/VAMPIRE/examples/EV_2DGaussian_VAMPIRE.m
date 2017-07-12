%==============================================================================
% This code is part of the VAMPIRE app for the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR/VAMPIRE.m 
%==============================================================================
%
% 2D Example for Multilevel Mass-Preserving Image Registration using VAMPIRE
% 
% 	- data                 synthetic 2D Gaussian blobs (level 3:8, full
%                          resolution: 256x256)
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             SSD
%   - pre-registration     none
%   - regularizer          mfHyperElastic
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

% finally: run the Mass-Preserving Non-Parametric Image Registration
NPIRpara            = optPara('NPIR-GN');
NPIRpara.lineSearch = @ArmijoDiffeomorphic;
NPIRpara.solver     = @VAMPIREsolveGN_PCG;

[yc,wc,his] = MLIR(ML, 'NPIRobj', @VAMPIRENPIRobjFctn, ...
    'parametric', false, 'NPIRpara', NPIRpara, 'minLevel', 5);

% [reg,regOptn] = regularizer('reset','regularizer','mfElastic',...
%   'alpha',alpha,'alphaLength',alphaLength,'alphaArea',alphaArea,...
%   'alphaVolume',alphaVolume);
% 
% % finally: run the  Non-Parametric Image Registration
% [yc,wc,his] = MLIR(ML, 'NPIRobj', @NPIRobjFctn, ...
%   'parametric', false, 'minLevel', 4);

%% Plot Results
% Compute resulting image: dataT(yc) * det(D(yc))
Topt = reshape(linearInter(dataT,omega,center(yc,m)) .* geometry(yc,m,'Jac','omega',omega), m);
figure;
subplot(2,2,1);
viewImage2D(dataT,omega,m,'colormap', 'gray(256)');
hold on; plotGrid(center(yc, m), omega, m, 'spacing', [5 5]); axis off; hold off;
title('Template (T) & Grid')
subplot(2,2,2);
viewImage2D(dataR,omega,m,'colormap', 'gray(256)'); axis off;
title('Reference (R)')
subplot(2,2,3);
viewImage2D(Topt,omega,m,'colormap', 'gray(256)'); axis off;
title('VAMPIRE result (Topt)')
subplot(2,2,4);
viewImage2D(abs(Topt-dataR),omega,m,'colormap', 'gray(256)'); axis off;
title('Absolute difference image of R and Topt');