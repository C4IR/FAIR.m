%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR:  Thin plate spline LM registration for hand example
%
% - load data (see setup2DhandData)
% - setup  viewer (viewImage2D), interpolator (splineInter), 
% - setup landmarks (LM)
% - run affine
%==============================================================================

clear, close all, help(mfilename)

setup2DhandData; 
xc = getCellCenteredGrid(omega,m);
cc = getTPScoefficients(LM(:,1:4),'theta',0);
[yc,yLM] = evalTPS(LM,cc,xc); 
LM(:,[5,6]) = reshape(yLM,[],2); 
P5_LM; % for nice plots
%==============================================================================
