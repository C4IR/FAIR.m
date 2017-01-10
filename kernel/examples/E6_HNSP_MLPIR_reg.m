%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: MLPIR, MultiLevel Parametric Image Registration, regularized
% 
%   - data                 HNSP, Omega=(0,2)x(0,1), level=3:7, m=[256,128]
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             SSD
%   - pre-registration     splineTransformation2D, p=[8,8], alpha=2e6, M=eye
%   - optimization         Gauss-Newton
% ===============================================================================

clear, close all, help(mfilename);

% load data, set viewer, interpolator, distance
setup2DHNSPData; 
imgModel('reset','imgModel','splineInter','regularizer','none','theta',0);
distance('reset','distance','SSD');

% setup transformation and regularization 
p = [8,8]; % using a p-spline grid
trafo('reset','trafo','splineTransformation2D','omega',omega,'m',m,'p',p);
w0 = trafo('w0'); wStop = w0;             % get starting point and stopping
M = 2e6*speye(length(w0)); wRef = w0;     % initilize the regularization

% run MLPIR
optn = {'M',M,'wRef',wRef,'Plots',@FAIRplots,'plotIter',0,'plotMLiter',0};
wc = MLPIR(ML,optn{:});
%==============================================================================
