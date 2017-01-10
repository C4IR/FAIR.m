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
%   - pre-registration     rigid2D
%   - optimization         Gauss-Newton
% ===============================================================================

clear, close all, help(mfilename);

setup2DHNSPData;                            % set up data
distance('reset','distance','SSD');       % specify distance measure
imgModel('reset','imgModel','splineInter');     % specify interpolator
trafo('reset','trafo','rigid2D');         % specify transformation
[wc,his] = MLPIR(ML,'minLevel',3,'plotMLiter',1);
%==============================================================================
