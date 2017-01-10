%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
% ##2
%==============================================================================
% 
%   - data                 3d-brains, Omega=(0,20)x(0,10)x(0,20), level=3:6, m=[128,64,128]
%   - viewer               imgmontage
%   - interpolation        splineInterMex
%   - distance             SSD
%   - pre-registration     affine3D
%   - regularizer          mfElastic
%   - optimizer            Gauss-Newton
% ===============================================================================

setup3DbrainData;

% configure the plotting routine
FAIRplots('clear')
Dshow = @(T,R,omega,m) viewIP(abs(T-R),omega,m,'colormap',gray(256));
FAIRplots('set','Dshow',Dshow);

% finally: run the MultiLevel Non-Parametric Image Registration, no pre-registration
[yc,wc,his] = MLIR(ML,'parametric',0);

%==============================================================================
