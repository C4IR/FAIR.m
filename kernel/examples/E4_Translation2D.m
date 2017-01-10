%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: translation of an US image
%
% - load data                  (setup2DUSData)
% - transform                  (translation2D)
% - interpolate                (linearInter) 
% - and visualize              (viewImage2D)
% 
%==============================================================================

setup2DUSData;  
fprintf('trafo=%s\n',translation2D([]));
wc = [-50;0]; 
xc = getCellCenteredGrid(omega,m);  
yc = translation2D(wc,xc);
Tc = linearInter(dataT,omega,yc);
FAIRfigure(2); viewImage2D(Tc,omega,m,'colormap','gray(256)'); 
%==============================================================================

