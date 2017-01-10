%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: 2D interpolation and transformations
%
% - load data                  (setup2DUSData)
% - generate holomorphic map
% - interpolate                (linearInter) 
% - and visualize              (viewImage2D)
% 
%==============================================================================

setup2DUSData; 
box = (omega(2:2:end)-omega(1:2:end));
xc = reshape(getCellCenteredGrid(omega,m),[],2);
yc = [(box(1)*((1 - 0.9*xc(:,2)/box(2)).*cos(pi*(1-xc(:,1)/box(1)))/2 + 0.5))
  (box(2)*(1-(1 - 0.9*xc(:,2)/box(2)).*sin(pi*(1-xc(:,1)/box(1)))))];

Tc = linearInter(dataT,omega,yc);
FAIRfigure(2); viewImage2D(Tc,omega,m,'colormap','gray(256)'); 