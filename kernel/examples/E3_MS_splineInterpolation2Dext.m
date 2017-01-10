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

setup2DUSData; m = 128*[3,2]; 
xc = getCellCenteredGrid(omega,m);
cc = @(reg,theta) getSplineCoefficients(dataT,'dim',2,'regularizer',reg,'theta',theta);

Tc   = @(theta) splineInter(cc('moments',theta),omega,xc);
name = @(theta) fullfile(FAIRpath,'temp',sprintf('US-MS-theta=%d.jpg',theta));

theta = 1e3
T = Tc(theta);
FAIRfigure(2); clf; 
viewImage2D(T,omega,m,'colormap','gray(256)');
imwrite(uint8(flipud(reshape(T,m)')),name(theta));
%==============================================================================
