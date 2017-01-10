%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: rotation of an US image
%
% - load data                  (setup2DUSData)
% - transform                  (rigid2D)
% - interpolate                (linearInter) 
% - and visualize              (viewImage2D)
%==============================================================================

setup2DUSData; 
fprintf('trafo=%s\n',rigid2D('[]'));
alpha = pi/6; R = [cos(alpha),-sin(alpha);sin(alpha),cos(alpha)];
center = (omega(2:2:end)+omega(1:2:end))'/2;
wc = [alpha;(eye(2)-R)*center]; 
xc = getCellCenteredGrid(omega,m);
yc = rigid2D(wc,xc);
Tc = linearInter(dataT,omega,yc);
FAIRfigure(2); viewImage2D(Tc,omega,m,'colormap','gray(256)'); 
%==============================================================================
