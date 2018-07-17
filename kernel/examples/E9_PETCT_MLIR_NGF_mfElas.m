%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
% ##2
%==============================================================================
% 
%   - data                 PETCT, Omega=(0,140)x(0,151), level=4:7, m=[128,128]
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             NGF
%   - pre-registration     affine2D
%   - regularizer          mfElastic
%   - optimization         Gauss-Newton
% ===============================================================================

clear; close all; help(mfilename)

setup2DPETCTData
theta = 0; % something to play with
imgModel('reset','imgModel','splineInter','regularizer','moments','theta',theta);
distance('reset','distance','NGF');
trafo('reset','trafo','affine2D'); wStop = trafo('w0'); w0 = wStop;
regularizer('reset','regularizer','mfElastic','alpha',0.05,'mu',1,'lambda',0);
regularizer('disp');

fprintf(' - run MLIR using sufficient amount of details (minLevel=4)\n');
yc =  MLIR(ML,'minLevel',4,'plotIter',0,'plotMLiter',0);

%==============================================================================
