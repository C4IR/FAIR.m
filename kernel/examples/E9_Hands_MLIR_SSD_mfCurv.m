%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% 
%   - data                 Hand, Omega=(0,20)x(0,25), 
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             SSD
%   - regularizer          mfCurvature
%   - optimizer            Gauss-Newton
% ===============================================================================

clear; close all; help(mfilename)

setup2DhandData
imgModel('reset','imgModel','splineInter','regularizer','moments','theta',0.01);
distance('reset','distance','SSD');
trafo('reset','trafo','affine2D');
regularizer('reset','regularizer','mfCurvature','alpha',1000);

yc = MLIR(ML,'parametric',1);
showResults(ML,yc);

%==============================================================================
