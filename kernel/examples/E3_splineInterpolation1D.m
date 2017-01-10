%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: spline interpolation in 1D, from the book
%
%==============================================================================

dataT = [0,2,2,2,1]; 
m     = length(dataT);  
omega = [0,m]; 
xc    = getCellCenteredGrid(omega,m);
B     = spdiags(ones(m,1)*[1,4,1],[-1:1],m,m);
T     = B\reshape(dataT,m,1);
xf    = linspace(-1,6,101);             
Tf    = splineInter(T,omega,xf);

FAIRfigure(1); clf; 
ph = plot(xc,dataT,'.',xf,Tf,'g-','markersize',30); 
%==============================================================================
