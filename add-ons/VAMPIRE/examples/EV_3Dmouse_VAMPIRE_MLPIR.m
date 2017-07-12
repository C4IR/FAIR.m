%==============================================================================
% This code is part of the VAMPIRE app for the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR/VAMPIRE.m 
%==============================================================================
%
% 3D Example for Parametric Multilevel Mass-Preserving Image Registration using VAMPIRE
% 
%   - data                 Cardiac gated PET images of a mouse heart.
%                          (level 4:6, full resolution: 40x40x40)
%                          Images show the heart in systole and diastole.
%   - viewer               viewImage3D
%   - interpolation        splineInterMex (regularizer='moments', theta=0.01)
%   - distance             SSD
%   - pre-registration     none
%   - transformation       splineTransformation3Dsparse
%   - regularizer          parametric FAIR regularizer
%   - optimizer            Gauss-Newton with ArmijoBacktrack linesearch
%
% Acknowledgement:
% ----------------
% Thanks to the European Institute for Molecular Imaging (EIMI) and SFB 656, 
% University of Muenster, Germany for supplying this interesting data.

setup3DmouseData

% prepare the plot
FAIRplots('clear')
Dshow = @(T,R,omega,m) viewIP(abs(T-R),omega,m,'colormap',gray(256));
FAIRplots('set','Dshow',Dshow);

% set the regularizer
alpha       = 1;
alphaLength = 1;
alphaArea   = 0.1;
alphaVolume = 2;
regularizer('reset', 'regularizer', 'mfHyperElastic', 'alpha', alpha, ...
    'alphaLength', alphaLength, ...
    'alphaArea',   alphaArea, ...
    'alphaVolume', alphaVolume);

% set the parametric transformation
p = [10 10 10];
trafo('reset','trafo','splineTransformation3Dsparse','omega',omega,'m',m,'p',p);




% finally: run the MultiLevel Non-Parametric Image Registration
PIRpara            = optPara('PIR-GN');
PIRpara.solver     = @VAMPIREsolveGN_PCG;

[wc,his] = MLPIR(ML, 'PIRobj', @VAMPIREPIRobjFctn, ...
    'getGrid', @getNodalGrid,'PIRpara',PIRpara);

%% Plot results
yc = trafo(wc,getNodalGrid(omega, m));

% Compute resulting image: dataT(yc) * det(D(yc))
Topt = reshape(linearInter(dataT,omega,center(yc,m)) .* geometry(yc,m,'Jac','omega',omega), m);

figure(1); clf;
subplot(2,2,1)
viewSlices(dataT,omega,m);
title('dataT, template');

subplot(2,2,2)
viewSlices(dataR,omega,m);
title('dataR, reference');

subplot(2,2,3)
viewSlices(Topt,omega,m);
title('dataT(y).*detDy');

subplot(2,2,4)
viewSlices((Topt-dataR),omega,m);
title('final residual');
