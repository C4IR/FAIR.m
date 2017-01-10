%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
% ##2
%==============================================================================
%
% - load data (see setup2DHNSPData)
% - setup  viewer          (viewImage2D), 
%           interpolator    (splineInter), 
%           distance        (SSD), 
%           transformation  (affine2D, not regularized)
%           regularizer     (curvature,matrix based)
% - run optimization 
% ===============================================================================

clear, close all, help(mfilename);

% load data, set viewer, interpolator, transformation, distance
setup2DHNSPData
imgModel('reset','imgModel','splineInter','regularizer','moments','theta',1e-2);
distance('reset','distance','SSD');
trafo('reset','trafo','affine2D');
regularizer('reset','regularizer','mbCurvature','alpha',1e1);

% run optimization
yc =  MLIR(ML,'maxLevel',7);

%==============================================================================
