%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
% ##2
%==============================================================================
% 
%  Example for 2D registration of hand data,
%
%   - data                 hands, omega=(0,20)x(0,25), level=3:7, m=[128,128]
%   - viewer               viewImage2D
%   - image model          splineInter
%   - distance             SSD
%   - pre-registration     affine2D
%   - regularizer          mfElastic
%   - optimization         Gauss-Newton
% ===============================================================================

close all, help(mfilename)

setup2DhandData
viewImage('reset','viewImage','viewImage2D','colormap',bone(256),'axis','off');
imgModel('reset','imgModel','splineInter','regularizer','moments','theta',1e-2);
trafo('reset','trafo','affine2D');
distance('reset','distance','SSD');
regularizer('reset','regularizer','mfElastic','alpha',1e3,'mu',1,'lambda',0);

[yc,wc,his] = MLIR(ML,'parametric',1,'plotIter',1,'plotMLiter',1);

showResults(ML,yc)

%==============================================================================
