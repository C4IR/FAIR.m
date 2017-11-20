%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: rotation of an US image, based on the book's version
%
% - load data                  (setup2DUSData)
% - transform                  (book_rigid2D)
% - interpolate                (linearInter) 
% - and visualize              (viewImage2D)
%==============================================================================

setup2DUSData; 

c   = (omega(2:2:end)+omega(1:2:end))'/2; alpha = pi/6; 
rot = [cos(alpha),-sin(alpha);sin(alpha),cos(alpha)];
wc  = [alpha;(eye(2)-rot)*c]; 
xc  = getCellCenteredGrid(omega,m);
yc  = rigid2D(wc,xc);
Tc  = linearInter(dataT,omega,yc);
FAIRfigure(2); viewImage2D(Tc,omega,m,'colormap','gray(256)'); 
%==============================================================================

 
