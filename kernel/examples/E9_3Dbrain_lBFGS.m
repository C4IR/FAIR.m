%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
% ##2
%==============================================================================
% 
%   - data                 3d-brains, Omega=(0,20)x(0,10)x(0,20), level=:6, m=[128,64,128]
%   - viewer               imgmontage
%   - interpolation        linearInter
%   - distance             SSD
%   - pre-registration     affine3D
%   - regularizer          mfElastic
%   - optimizer            lBFGS
% ===============================================================================

setup3DbrainData;

% prepare the plot
FAIRplots('clear')
Dshow = @(T,R,omega,m) viewIP(abs(T-R),omega,m,'colormap',gray(256));
FAIRplots('set','Dshow',Dshow);

% finally: 
% run the MultiLevel Non-Parametric Image Registration
% with lBFGS optimizer


NPIRpara  = optPara('lBFGS');
[yc,wc,his] = MLIR(ML,'NPIRobj',@NPIRBFGSobjFctn,'NPIRpara',NPIRpara);

%==============================================================================
