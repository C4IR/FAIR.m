%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
% ##2
%==============================================================================
%
% Matrix-free Hyperelastic registration of 3D brain data
% 
%   - data                 3D brain, Omega=(0,20)x(0,10)x(0,20), 
%                          level=3:5, m=[128,64,128]
%   - viewer               imgmontage
%   - image model          splineInterMex
%   - distance             SSD
%   - pre-registration     none
%   - regularizer          mfHyperElastic
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

setup3DbrainData;

% prepare the plot
FAIRplots('clear')
Dshow = @(T,R,omega,m) viewIP(abs(T-R),omega,m,'colormap',gray(256));
FAIRplots('set','Dshow',Dshow);

regularizer('reset','regularizer','mfHyperElastic',...
  'alpha',100,...
  'alphaLength',1,...
  'alphaArea',0.1,...
  'alphaVolume',1);

%% finally: run the MultiLevel Non-Parametric Image Registration
NPIRpara    = optPara('NPIR-GN');
NPIRpara.lineSearch = @ArmijoDiffeomorphic;

[yc,wc,his] = MLIR(ML,'parametric',0,'NPIRpara',NPIRpara,'maxLevel',5);
%==============================================================================

