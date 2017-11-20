%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
% ##2
%==============================================================================
% 
%   - data                 Hand, Omega=(0,20)x(0,25), 
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             SSD
%   - transformation       affine2D
%   - regularizer          mfElastic
%   - optimizer            Gauss-Newton
% ===============================================================================

close all, help(mfilename);

% load data, set viewer, interpolator, transformation, distance
setup2DHNSPData
imgModel('reset','imgModel','splineInter','regularizer','moments','theta',1e-2);
distance('reset','distance','SSD');
trafo('reset','trafo','affine2D');
regularizer('reset','regularizer','mfElastic','alpha',5e2,'mu',1,'lambda',0);

% run MLIR
[yc,wc,his] = MLIR(ML,'maxLevel',8);

%==============================================================================
