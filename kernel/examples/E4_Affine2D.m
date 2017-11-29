%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: 2D interpolation and transformations
%
% - load data                  (setup2DUSData)
% - generate affine linear map (affine2D)
% - interpolate                (linearInter) 
% - and visualize              (viewImage2D)
% 
%==============================================================================

setup2DUSData;  
wc = [1 -0.2 50, 0, 0.75 50]'; 
xc = getCellCenteredGrid(omega,m);  
yc = affine2D(wc,xc);
Tx = linearInter(dataT,omega,xc);
Ty = linearInter(dataT,omega,yc);
FAIRfigure(1); viewImage2D(Tx,omega,m,'colormap','gray(256)'); 
FAIRfigure(2); viewImage2D(Ty,omega,m,'colormap','gray(256)'); 
%==============================================================================
