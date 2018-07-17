%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
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
%   - pre-registration     rotation2D
%   - optimization         Gauss-Newton
% ===============================================================================

clear, close all, help(mfilename);

setup2DHNSPData;                            % set up data
distance('reset','distance','SSD');       % specify distance measure
imgModel('reset','imgModel','splineInter');   % specify interpolator and
center = (omega(2:2:end)-omega(1:2:end))'/2;
trafo('reset','trafo','rotation2D','c',center); 
[wc,his] = MLPIR(ML,'minLevel',3,'plotMLiter',1);
%==============================================================================
