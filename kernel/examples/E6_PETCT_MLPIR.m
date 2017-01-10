%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Example for MLPIR, MultiLevel Parametric Image Registration
% 
%   - data                 PETCT, Omega=(0,140)x(0,151), level=4:7, m=[128,128]
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             MI
%   - pre-registration     affine2D
%   - optimization         Gauss-Newton
% ===============================================================================

clear, close all, help(mfilename);

% load data, set viewer
setup2DPETCTData;

% setup interpolator, distance
theta = 1e2; % something to play with
imgModel('reset','imgModel','splineInter','regularizer','moments','theta',theta);
distance('reset','distance','MI'); 

% setup transformation, regularization 
trafo('reset','trafo','affine2D'); 
wc =  MLPIR(ML,'plotIter',0,'plotMLiter',0);
%==============================================================================
