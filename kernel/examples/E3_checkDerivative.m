%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: derivative check (use data from E3_splineInterpolation2D)
%
%==============================================================================


% setup test data
dataT = flipud([1,2,3,4;1,2,3,4;4,4,4,4])'; 
m     = size(dataT); 
omega = [0,m(1),0,m(2)]; 
B     = @(i) spdiags(ones(m(i),1)*[1,4,1],[-1:1],m(i),m(i));
T     = B(1)\dataT/B(2);
xf    = reshape(getCellCenteredGrid(omega,10*m),[],2);

fctn = @(x) splineInter(T,omega,x);
figure(1); clf;
[fig,ph,th] = checkDerivative(fctn,xf(:),'fig',1);
%==============================================================================
