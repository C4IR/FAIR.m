%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: 2D imgModelpolation and visualization
%
% the tutorial explores the scale-space idea for 2D data
% the spline based model is based on the minimization of
% D(T)+theta S(T)!= min, 
% where D(T) is a data fitting term and S(T) is the linearized bending energy
%
% - load data ('US.jpg')
% - display and visualize data  (viewImage2D)
% - compute approximations (spline, various theta's) and visualize
% 
%==============================================================================

clear, close all, help(mfilename); echo on

%% load data, define a doman and an initial discretization
dataT = double(imread('US.jpg'));
omega = [0,size(dataT,1),0,size(dataT,2)];
m     = [128,128]/2;
xc    = getCellCenteredGrid(omega,m);

% setup image viewer
viewImage('reset','viewImage','viewImage2D','colormap','gray(256)');

% setup spline imgModelpolation
imgModel('reset','imgModel','splineInter','regularizer','moments');
imgModel('disp');

scaleParameter = [logspace(3,-2,23),0];
titleStr = @(j) title(...
  sprintf('scale-space, \\theta=%s',num2str(scaleParameter(j))),'fontsize',30);

for j=1:length(scaleParameter);
  % set scale parameter and compute spline coefficients
  imgModel('set','theta',scaleParameter(j));
  T = imgModel('coefficients',dataT,[],omega,'out',0);
  
  % imgModelpolate image
  Tc = imgModel(T,omega,xc);
  
  if j==1, % initilize figure
    figure(1); clf;
    vh = viewImage(Tc,omega,m); titleStr(j);
    FAIRpause;
  else
    set(vh,'cdata',reshape(Tc,m)'); titleStr(j);
    FAIRpause(1/500);
  end;
end;
  
echo off
%==============================================================================
