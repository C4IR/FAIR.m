%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: spline interpolation in 2D, from the book
%
%==============================================================================

dataT = flipud([1,2,3,4;1,2,3,4;4,4,4,4])'; 
m     = size(dataT); 
omega = [0,m(1),0,m(2)]; 
M     = {m,10*m};   % two resolutions
xc    = reshape(getCellCenteredGrid(omega,M{1}),[],2);
xf    = reshape(getCellCenteredGrid(omega,M{2}),[M{2},2]);

B  = @(i) spdiags(ones(m(i),1)*[1,4,1],[-1:1],m(i),m(i));
T  = B(1)\dataT/B(2);
Tc = splineInter(T,omega,xc(:));
Tf = splineInter(T,omega,xf(:));

FAIRfigure(1); clf;
ph  = plot3(xc(:,1),xc(:,2),Tc(:),'ro');    hold on;
qh = surf(xf(:,:,1),xf(:,:,2),reshape(Tf,M{2}));
%==============================================================================
