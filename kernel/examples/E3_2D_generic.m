%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: 2D interpolation and visualization
%
% generic example for usage of interpolation and visualization
% loads data, initializes viewer and interpolator, shows some results
%
% - load data ('US.jpg')
% - display and visualize data  (viewImage2D)
% - view image in high resolution and low resolution
%==============================================================================

clear, close all, help(mfilename); 

fprintf('%s\n','generic interpolation and visualization')

fprintf('%s\n','load data, setup image viewer')
dataT = double(imread('US.jpg'));
m     = size(dataT);
omega = [0,size(dataT,1),0,size(dataT,2)];

% setup image viewer
viewImage('reset','viewImage','viewImage2D','colormap',gray(256));

% shortcuts for labeling off plots
dimstr  = @(m) sprintf('%s=[%s]',inputname(1),sprintf(' %d',m));
lblstr  = @(str,m) sprintf('%s, %s',str,dimstr(m));

fprintf('%s\n','visualize original data')
FAIRfigure(1,'figname',mfilename); clf; 
viewImage(dataT,omega,m);
title(lblstr('highres',m),'fontsize',20); set(gca,'fontsize',20); FAIRpause;

fprintf('%s\n','setup interpolator')
imgModel('reset','imgModel','splineInter','regularizer','moments','theta',1e-2);
T = getSplineCoefficients(dataT,'out',0);

fprintf('%s\n','generate some points and interpolate')
m  = [256,128]/4;
xc = getCellCenteredGrid(omega,m);
Tc = imgModel(T,omega,xc);

FAIRfigure(2,'figname',mfilename); clf; 
viewImage(Tc,omega,m);
title(lblstr('lowres',m),'fontsize',20);  set(gca,'fontsize',20);

%==============================================================================
