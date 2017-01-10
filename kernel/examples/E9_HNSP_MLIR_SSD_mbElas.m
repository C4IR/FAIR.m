%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
% ##2
%==============================================================================
% 
%   - data                 HNSP, Omega=(0,2)x(0,1), m=[ 512 256]
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             SSD
%   - pre-registration     affine2D
%   - regularizer          mbElastic
% ===============================================================================

close all, help(mfilename);

setup2DHNSPData
imgModel('reset','imgModel','splineInter','regularizer','moments','theta',1e-2);
distance('reset','distance','SSD');
trafo('reset','trafo','affine2D');
regularizer('reset','regularizer','mbElastic','alpha',5e2,'mu',1,'lambda',0);
% regularizer('reset','regularizer','mfElastic','alpha',5e2,'mu',1,'lambda',0);
% regularizer('reset','regularizer','mbCurvature','alpha',1e1);
% regularizer('reset','regularizer','mfCurvature','alpha',1e1);

[yc,wc,his] = MLIR(ML,'maxLevel',7,'parametric',1);

%==============================================================================
