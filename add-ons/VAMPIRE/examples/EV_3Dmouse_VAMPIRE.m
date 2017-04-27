% 3D Example for Multilevel Mass-Preserving Image Registration using VAMPIRE
% 
% (c) Fabian Gigengack and Lars Ruthotto 2011/02/04, see FAIR.2 and FAIRcopyright.m.
% http://www.uni-muenster.de/EIMI/
% http://www.mic.uni-luebeck.de/
%
%   - data                 Cardiac gated PET images of a mouse heart.
%                          (level 4:6, full resolution: 40x40x40)
%                          Images show the heart in systole and diastole.
%   - viewer               viewImage3D
%   - interpolation        splineInterMex (regularizer='moments', theta=0.01)
%   - distance             SSD
%   - pre-registration     none
%   - regularizer          mfHyperElastic
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

% finally: run the MultiLevel Non-Parametric Image Registration
NPIRpara            = optPara('NPIR-GN');
NPIRpara.lineSearch = @ArmijoDiffeomorphic;
NPIRpara.solver     = @VAMPIREsolveGN_PCG;

[yc,wc,his] = MLIR(ML, 'NPIRobj', @VAMPIRENPIRobjFctn, ...
    'parametric', false, 'NPIRpara', NPIRpara);

%% Plot results
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
