%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: multiscale spline interpolation in 2D
%
% - load data   (setup2DUSData)
% - interpolate (splineInter) enabling different theta's
% - visualize   (viewImage2D)
% 
%==============================================================================

setup2DUSData; m = 128*[3,2]; xc = getCellCenteredGrid(omega,m);
T = getSplineCoefficients(dataT,'dim',2,'regularizer','gradient','theta',50);
figure(2); clf; viewImage2D(splineInter(T,omega,xc),omega,m); 
colormap(gray(256));
%==============================================================================
