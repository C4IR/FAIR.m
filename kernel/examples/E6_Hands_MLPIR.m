%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: MLPIR, MultiLevel Parametric Image Registration
% 
%   - data                 Hand, Omega=(0,20)x(0,25), level=3:7, m=[128,128] 
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             SSD
%   - pre-registration     affine2D, not regularized
%   - optimization         Gauss-Newton
% ===============================================================================

clear, close all, help(mfilename);

% load data, initialize image viewer, interpolator, transformation, distance
setup2DhandData
imgModel('reset','imgModel','splineInter','regularizer','moments','theta',1e-1);
distance('reset','distance','SSD');
trafo('reset','trafo','affine2D');

% run Multilevel Parametric Image Registration 
wc = MLPIR(ML,'plotIter',0,'plotMLiter',1);
%==============================================================================
