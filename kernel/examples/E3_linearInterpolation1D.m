%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: linear interpolation in 1D, from the book
%
%==============================================================================

TD    = [1,2,3,3]; 
m     = length(TD); 
omega = [0,4]; TD = reshape(TD,m,1);
xc    = getCellCenteredGrid(omega,m);
Tc    = linearInter(TD,omega,xc);
xf    = linspace(omega(1)-1,omega(2)+1,201);
Tf    = linearInter(TD,omega,xf);
FAIRfigure(1); clf; 
ph = plot(xc,TD,'b+',xc,Tc,'ro',xf,Tf,'k-');
%==============================================================================
