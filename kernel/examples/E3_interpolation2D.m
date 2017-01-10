%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: interpolation in 2D
%
% - load data   (setup2DUSData)
% - interpolate (linearInter) on different resolutions
% - visualize   (viewImage2D)
% - image model (splineInter) on different resolutions
% - visualize   (viewImage2D)
%==============================================================================

setup2DUSData; close all; 

T = dataT; 
xc = @(m) getCellCenteredGrid(omega,m); 

imgModel('set','imgModel','linearInter');
for p=5:7,
  m = 2^p*[1,1]; 
  Tc = imgModel(T,omega,xc(m));
  FAIRfigure(p-4); viewImage2D(Tc,omega,m); colormap(gray(256));
end;

imgModel('set','imgModel','splineInter');
T = getSplineCoefficients(dataT,'regularizer','moments','theta',100);
for p=5:7,
  m = 2^p*[1,1]; 
  Tc = imgModel(T,omega,xc(m));
  figure(p-1); viewImage2D(Tc,omega,m); colormap(gray(256));
end;
%==============================================================================
