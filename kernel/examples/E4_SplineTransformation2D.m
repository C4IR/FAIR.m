%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: spline transformation of an US image
%
% - load data                  (setup2DUSData)
% - transform                  (splineTransformation2D)
% - interpolate                (linearInter) 
% - and visualize              (viewImage2D)
% 
%==============================================================================

setup2DUSData; 
p  = [5,4]; 
xc = getCellCenteredGrid(omega,m);  
splineTransformation2D([],xc,'omega',omega,'m',m,'p',p);  
w1 = zeros(p); w2 = zeros(p);  w2(3,2) = 3; 
wc = [w1(:);w2(:)]; 
yc = splineTransformation2D(wc,xc);
Tc = linearInter(dataT,omega,yc);
FAIRfigure(2); viewImage2D(Tc,omega,m,'colormap','gray(256)'); 
%==============================================================================

