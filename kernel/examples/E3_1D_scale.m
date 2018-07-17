%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: 1D image model, interpolation, scale-space in approximations
%
% the tutorial explores the scale-space idea for 1D data
% the spline based model is based on the minimization of
% D(T)+theta S(T)!= min, 
% where D(T) is a data fitting term and S(T) is the linearized bending energy
%
% - setup 1D data
% - visualize data
% - compute approximations (spline, various theta's) and visualize
%==============================================================================

clear, close all, help(mfilename); 
setBreakpoint(mfilename,{27,39,69});

%% setup data 
fprintf('%s\n','generate noisy data ');

omega = [0,10]; m = 21;
dataX = getCellCenteredGrid(omega,m);
dataT  = rand(m,1); dataT([1,end]) = 0; % T should be compactly supported

% visualize the data
FAIRfigure(2,'position',[2620 867 560 420]); clf;
ph = plot(dataX,0*dataT,'k.',dataX,dataT,'b.','markersize',30); hold on; 
title('1D multiscale','fontsize',20);
axis([omega(1)-1,omega(2)+1,min(dataT)-1,max(dataT)+1]);
lstr = {'location','data'}; legend(ph,lstr,'location','northwest');

%% discretization for the model  -----------------------------------------------
m  = 101; % discretization for the model, h=omega./m
xc = getCellCenteredGrid(omega,m);

% the coefficients depend on the scale-parameter theta
T = @(theta) getSplineCoefficients(dataT,...
  'regularizer','moments','theta',theta,'dim',1,'out',0);

% the continuous model
Tc = @(theta) splineInter(T(theta),omega,xc);

% show results for various theta's: fine to coarse scale
fprintf('%s\n','show results for various theta''s');
theta = [0,logspace(-3,3,51)];
for j=1:length(theta),
  fprintf('.');
  if j>1, 
    set(qh,'visible','off');  
  end;
  qh = plot(xc,Tc(theta(j)),'g-','linewidth',3);
  ylabel(sprintf('\\theta=%s',num2str(theta(j))),'fontsize',20); 
  FAIRpause(1/10)
  if j==1,  
    lstr{3} = 'spline';  
    legend([ph;qh],lstr,'location','northwest'); 
    FAIRpause; 
  end;
end;
fprintf('\n'); 

% show results for various theta's: coarse to fine scale
for j=length(theta):-1:1,
  fprintf('.');
  set(qh,'visible','off');
  qh = plot(xc,Tc(theta(j)),'g-','linewidth',3);
  ylabel(sprintf('\\theta=%s',num2str(theta(j))),'fontsize',20); 
  FAIRpause(1/10)
  if j==1, FAIRpause; end;
end;
fprintf('\n'); 

%==============================================================================
