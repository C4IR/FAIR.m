%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
% ##2
%==============================================================================
% 
%   - data                 PETCT, Omega=(0,140)x(0,151), level=4:7, m=[128,128]
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             MI
%   - pre-registration     rigid2D
%   - regularizer          mbCurvature
%   - optimization         lBFGS
% ===============================================================================

close all, help(mfilename);

setup2DPETCTData
imgModel('reset','imgModel','splineInter','regularizer','none','theta',1e-3);
distance('reset','distance','MI','nT',32,'nR',32);
trafo('reset','trafo','rigid2D');
regularizer('reset','regularizer','mfCurvature','alpha',1e-1);

PIRpara = optPara('lBFGS','solver','backslash');
NPIRpara = optPara('lBFGS','solver',regularizer('get','solver'),'maxIter',40);

[yc,wc,his] = MLIR(ML,'PIRobj',@PIRBFGSobjFctn,'PIRpara',PIRpara,...
  'NPIRobj',@NPIRBFGSobjFctn,'NPIRpara',NPIRpara,...
  'minLevel',4,'maxLevel',7,'parametric',1,'plotMLiter',0);

%==============================================================================
