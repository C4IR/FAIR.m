%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
% ##2
%==============================================================================
%
% Hyperelastic registration of 2D images disc "o" to image "C"
% 
%   - data                 disc" and "C", Omega=(0,1)x(0,1), level=3:7, m=[128,128]
%   - viewer               viewImage2D
%   - image model          splineInter
%   - distance             SSD
%   - pre-registration     none
%   - regularizer          mbHyperElastic
%   - optimizer            Gauss-Newton with ArmijoDiffeomorphic 
%
% For more info see:
%
% @article{BurgerEtAl2013,
%    author = {Burger, M and Modersitzki, J and Ruthotto, L},
%    title = {{A hyperelastic regularization energy for image registration}},
%    journal = {SIAM Journal on Scientific Computing},
%    year = {2013},
%    volume = {35},
%    number = {1},
%    pages = {B132--B148},
%    keywords = {Image Registration},
%    doi = {10.1137/110835955},
% }
%==============================================================================

clear, close all, help(mfilename), 
setup2Ddisc2CData;

% prepare the plot
FAIRplots('clear')
% Dshow = @(T,R,omega,m) viewIP(abs(T-R),omega,m,'colormap',gray(256));
% FAIRplots('set','Dshow',Dshow);

% initialize the regularizer for the non-parametric part
alpha       = 1;
alphaLength = 100;
alphaArea   = 0;
alphaVolume = 18;

[reg,regPara] = regularizer('reset','regularizer','mbHyperElastic',...
  'alpha',alpha,...
  'alphaLength', alphaLength,...
  'alphaArea',alphaArea,...
  'alphaVolume',alphaVolume);

%% finally: run the MultiLevel Non-Parametric Image Registration
NPIRpara            = optPara('NPIR-GN');
NPIRpara.lineSearch = @ArmijoDiffeomorphic;
NPIRpara.maxIter    = 50;

[yc,wc,his] = MLIR(ML, 'parametric', false,...
    'minLevel', 3, 'maxLevel', 7, 'NPIRpara', NPIRpara,'plots',1);
%==============================================================================
	