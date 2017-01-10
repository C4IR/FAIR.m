%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: using MATLAB's linear interpolation in 1D 
%
%==============================================================================

T     = [1,2,3,3]'; 
m     = length(T);
omega = [0,4]; 

xc = getCellCenteredGrid(omega,m);
Tc = linearInterMatlab(T,omega,xc);
xf = linspace(omega(1)-1,omega(2)+1,201);
Tf = linearInterMatlab(T,omega,xf);
clf; ph = plot(xc,0*Tc,'b+',xc,Tc,'ro',xf,Tf,'k-');

fctn = @(x) linearInterMatlab(T,omega,x);
fig = checkDerivative(fctn,xc+1e-2); FAIRpause(2); close(fig);
%==============================================================================
