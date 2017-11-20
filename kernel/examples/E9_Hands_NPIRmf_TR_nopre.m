%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
% ##2
%==============================================================================
%
% - initialize
%     data:            Hands,           see setup2DhandData, level = 5;
%     visualization:   viewImage2D,     see viewImage.m
%     interpolation:   splineInter,   see inter.m
%     distance:        SSD,             see distance.m
%     regularizer:     mfElastic,       see regularizer.m
% - initialize FAIR plots
% - run optimization
%     NPIR:            Non-Parametric Image Registration, Trust-Region
% ===============================================================================

clear; close all; help(mfilename)

setup2DhandData
level = 5; omega = ML{level}.omega; m = ML{level}.m;
imgModel('reset','imgModel','splineInter','regularizer','moments','theta',0.01);
distance('reset','distance','SSD');
regularizer('reset','regularizer','mfElastic','alpha',1000,'mu',1,'lambda',0);

[T,R] = imgModel('coefficients',ML{level}.T,ML{level}.R,omega);
xc    = getStaggeredGrid(omega,m); % starting guess and reference for regularization
Rc    = imgModel(R,omega,center(xc,m)); 

% - initialize FAIR plots
FAIRplots('set','mode','TR-mf','fig',1);
FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m)); 

% - run Non-Parametric Image Registration (Gauss-Newton)
fctn = @(yc) NPIRobjFctn(T,Rc,omega,m,xc,yc);  fctn([]);
yc = TrustRegion(fctn,xc,'Plots',@FAIRplots);
%==============================================================================
