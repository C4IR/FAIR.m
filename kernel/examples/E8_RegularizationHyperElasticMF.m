%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for Regularization: Check Derivative of mfHyperElastic
%
% where
%  alpha regularization parameter, weights regularization versus 
%        distance in the joint objective function, alpha = 1 here
%  yRef  is a reference configuration, e.g. yRef = x or a 
%        pre-registration result
%  B     a discrete partial differential operator either in explicit
%        matrix form or as a structure containing the necessary
%        parameters to compute B*y
% see also regularizer E8_regularization_MB
%==============================================================================

clear, close all, help(mfilename);


% initialize the regularization and create a starting  point
alphaLength = 1;
alphaArea   = 1;
alphaVolume = 1;
regularizer('reset','regularizer','mfHyperElastic','alpha',1,'alphaLength',1,...
    'alphaArea',alphaArea,'alphaVolume',alphaVolume);
y0 = @(omega,m) randn((length(omega)/2)*prod(m+1),1);

% 2D example, initialize physical domain and number of discretization points
omega = [0,1,0,1]; m = [16,12];   % 

% test derivative of 2D implementation
fctn = @(yc) regularizer(yc,omega,m);  
checkDerivative(fctn,y0(omega,m));

% 3D example, initialize physical domain and number of discretization points
omega = [0,1,0,1,0,1]; m  = [16,12,8];

% test derivative of 3D implementation
fctn = @(yc) regularizer(yc,omega,m); 
checkDerivative(fctn,y0(omega,m));
%==============================================================================
