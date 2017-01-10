%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: 2D interpolation, handling of derivatives
%
% the tutorial explores the derivative of the interpolant 
% for 'linearInterMatlab', 'linearInter', 'splineInter'
%
% - setup data
% - display and visualize data
% - compute interpolants and check derivative
%==============================================================================

clear, close all, help(mfilename); echo on

%% exploring the derivatives of different interpolants for data given by USfair.jpg
fprintf('%s\n','load data')
dataT = double(imread('US.jpg'));
m     = floor(size(dataT)/4);
omega = [0,size(dataT,1),0,size(dataT,2)];
xc    = getCellCenteredGrid(omega,m);

%% these are the methods to be compared
methods = {'linearInterMatlab','linearInter','splineInter'};
for j=1:length(methods),
  
  % reset the interpolation scheme and compute coefficients (if necessary)
  imgModel('reset','imgModel',methods{j});
  T  = imgModel('coefficients',dataT,[],omega,'out',0);
  Tc = @(x) imgModel(T,omega,x);

  fprintf('testing %s, using regular x\n',imgModel)
  [f,ph,th] = checkDerivative(Tc,xc);

  % make plots nice
  set(gca,'fontsize',20,'YAxisLocation','right');
  set(ph,'linewidth',2); set(th,'fontsize',20);
  ylabel(sprintf('[%s], based on x',imgModel))
  FAIRpause;

  fprintf('test %s, using X+randn(size(X))\n',imgModel)
  [f,ph,th] = checkDerivative(Tc,xc+randn(size(xc)));

  % make plots nice
  set(gca,'fontsize',20,'YAxisLocation','right');
  set(ph,'linewidth',2); set(th,'fontsize',20);
  ylabel(sprintf('[%s], based on X+noise',imgModel));
  FAIRpause;
end;
echo off
%==============================================================================
