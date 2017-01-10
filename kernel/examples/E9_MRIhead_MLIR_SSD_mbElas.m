%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
% ##2
%==============================================================================
% 
%   - data                 MRI (head), Omega=(0,128)x(0,128), level=4:7, m=[128,128]
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             SSD
%   - pre-registration     rigid2D
%   - regularizer          mbElastic
%   - optimization         Gauss-Newton
% version 2015/05/20
% ===============================================================================


close all, help(mfilename)

setup2DMRIData
theta = 1e-1; % something to play with
imgModel('reset','imgModel','splineInter','regularizer','moments','theta',theta);
distance('reset','distance','SSD'); 
trafo('reset','trafo','rigid2D');
wStop = trafo('w0'); w0 = wStop;
regularizer('reset','regularizer','mbElastic','alpha',1e4,'mu',1','lambda',0');
regularizer('disp');

fprintf(' - run MLIR using sufficient amount of details (minLevel=4)\n');
yc =  MLIR(ML,'minLevel',4,'maxLevel',7,'plotIter',0,'plotMLiter',0);

%==============================================================================

