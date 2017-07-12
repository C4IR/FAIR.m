% 2D Example for Multilevel Mass-Preserving Image Registration using VAMPIRE
% 
% (c) Fabian Gigengack and Lars Ruthotto 2011/02/04, see FAIR.2 and FAIRcopyright.m.
% http://www.uni-muenster.de/EIMI/
% http://www.mic.uni-luebeck.de/
% 
%   - data                 synthetic 2D Gaussian blobs (level 3:8, full
%                          resolution: 256x256), data disturbed with
%                          Gaussian noise
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             SSD
%   - pre-registration     none
%   - regularizer          mfHyperElastic
%   - optimizer            Gauss-Newton with ArmijoBacktrack linesearch

setup2DGaussianData;

nT    = randn(size(dataT));
dataT = dataT + .4*norm(dataT(:))*nT/norm(nT(:));
nR    = randn(size(dataR));
dataR = dataR + .4*norm(dataR(:))*nR/norm(nR(:));

ML    = getMultilevel({dataT,dataR},omega,m);
%%
% prepare the plot
FAIRplots('clear')
Dshow = @(T,R,omega,m) viewIP(abs(T-R),omega,m,'colormap',gray(256));
FAIRplots('set','Dshow',Dshow);

% initialize the regularizer for the non-parametric part
alpha       = 100;
alphaLength = 1;
alphaArea   = 0;
alphaVolume = 1;
[reg,regOptn] = regularizer('reset','regularizer','mfHyperElastic',...
  'alpha',alpha,'alphaLength',alphaLength,'alphaArea',alphaArea,...
  'alphaVolume',alphaVolume);

% finally: run the MultiLevel Non-Parametric Image Registration
NPIRpara            = optPara('NPIR-GN');
NPIRpara.lineSearch = @ArmijoDiffeomorphic;
NPIRpara.solver     = @VAMPIREsolveGN_PCG;

[yc,wc,his] = MLIR(ML, 'NPIRobj', @VAMPIRENPIRobjFctn, ...
  'parametric', false, 'minLevel', 5, 'NPIRpara', NPIRpara);

%% Plot results
% Compute resulting image: dataT(yc) * det(D(yc))
Topt = reshape(linearInter(dataT,omega,center(yc,m)) .* ...
    geometry(yc,m,'Jac','omega',omega), m);
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
title('Absolute difference image of R and Topt')