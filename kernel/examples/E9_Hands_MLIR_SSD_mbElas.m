%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% 
%   - data                 Hand, Omega=(0,20)x(0,25), 
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             SSD
%   - regularizer          mbElastic
%   - optimizer            Gauss-Newton
% ===============================================================================

close all, help(mfilename)

% load data, initialize image viewer, interpolator, transformation, distance
setup2DhandData
viewImage('reset','viewImage','viewImage2D','colormap',bone(256),'axis','off');
imgModel('reset','imgModel','splineInter','regularizer','moments','theta',1e-2);
trafo('reset','trafo','affine2D');
distance('reset','distance','SSD');
regularizer('reset','regularizer','mbElastic','alpha',1e3,'mu',1,'lambda',0);


PIRobj = @PIRobjFctn;
PIRopt = optPara('GN');

NPIRobj = @NPIRobjFctn;
NPIRopt = optPara('NPIR-GN')

[yc,wc,his] = MLIR(ML,'PIRobj',PIRobj,'PIRopt',PIRopt,'NPIRobj',NPIRobj,...
  'NPIRopt',NPIRopt,'pause',1,'parametric',1);
  
% ===============================================================================
