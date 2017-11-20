%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
% ##2
%==============================================================================
% 
%   - data                 HNSP, Omega=(0,2)x(0,1), level=3:6, m=[512,256]
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             SSD
%   - pre-registration     affine2D
%   - regularizer          mbElastic
%   - optimization         TrustRegion
% ===============================================================================

clear; close all; help(mfilename)

setup2DHNSPData
imgModel('reset','imgModel','splineInter','regularizer','moments','theta',0.01);
distance('reset','distance','SSD');
trafo('reset','trafo','affine2D');
regularizer('reset','regularizer','mfElastic','alpha',1000,'mu',1,'lambda',0);

NPIRpara = optPara('TrustRegion')
yc = MLIR(ML,'maxLevel',6,'NPIRpara',NPIRpara);

%==============================================================================
