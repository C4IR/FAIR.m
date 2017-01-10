%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: 2D interpolation, basics ideas
%
% basic 2D interpolation example, initializing small data, plotting the data and a 
% visualization of the interpolant
%
% - setup data
% - display and visualize data
% - interpolate the data and visualize 
%==============================================================================

clear, close all, help(mfilename); echo on

% setup data 
fprintf('%s\n','generate data ...');

omega = [0,8,0,6];  % domain omega = (omega(1),omega(2)) x (omega(3),omega(4))
dataT = [           % the data given as m(1)-by-m(2) array
  4     1     1
  4     2     2
  4     3     3
  4     4     4 ];
m  = size(dataT);   % data size
xc = getCellCenteredGrid(omega,m);

% create shortcut for i-th components of the grid points
dim = size(omega,2)/2;
e = @(i) ((1:dim)' == i);  % i-th unit vector
x = @(i) reshape(xc,[],2)*e(i);

fprintf('[x1;x2;T]=\n');
disp(reshape([x(1);x(2);dataT(:)],[],3)');
FAIRpause;

% plot the data
fprintf('%s\n','visualize the data')
FAIRfigure(1); clf;
ph = plot3(x(1),x(2),dataT(:),'ro');
set(ph,'markersize',20,'markerfacecolor','r' ); hold on;
axis([omega,-1,5])
ah = [xlabel('$x^1$'),ylabel('$x^2$'),zlabel('$x^3$')];
set(ah,'interpreter','latex','fontsize',40);
set(gca,'fontsize',30); 

% interpolate the data and visualize 
mfine = 10*m;
Xfine = reshape(getCellCenteredGrid(omega,mfine),[mfine,2]);
T = getSplineCoefficients(dataT,'out',0); % spline coefficients
Tfine = splineInter(T,omega,Xfine(:));  % sample interpolant

sh = surf(Xfine(:,:,1),Xfine(:,:,2),reshape(Tfine,mfine));
set(sh,'facealpha',0.5)
FAIRpause; 

% visualize
fprintf('%s\n','alternative visualization using viewImage2D')
FAIRfigure(2); clf;
subplot(1,2,1); viewImage2D(dataT,omega,m,'colormap',gray(6)); hold on;
plotGrid(getNodalGrid(omega,m),omega,m)
title('visualization of data','fontsize',30)
set(gca,'fontsize',30)
subplot(1,2,2); viewImage2D(Tfine,omega,mfine,'colormap',gray(4)); hold on;
title('visualization of interpolant','fontsize',30)
set(gca,'fontsize',30)

echo off
%==============================================================================
