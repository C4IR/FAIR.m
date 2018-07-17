%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% 
% Tutorial for FAIR: shows how to visulaize 3D data with FAIR.
%
% - load data             (here 3D MRI)
% - setup viewer          (viewImage2D)
% - view  data
%==============================================================================

clear, %close all, 
help(mfilename), 

% setup data 
fprintf('%s\n','generate data (here: MRI)');
load mri;
omega = [0,20,0,20,0,10];    % specify physical domain
T     = double(squeeze(D));  
m     = size(T);             % get size of data

%%
fprintf('%s\n','use image montage for visualization of data')
FAIRfigure(1); clf; 
imgmontage(T,omega,m,'colormap',gray(100),'numbering','on');
title('visualization of 3D data - imgmontage','fontsize',30)

%%
fprintf('%s\n','use a volumetric view')
FAIRfigure(2); clf; 
volView(T,omega,m,'facecolor',[240,140,100]/256,'facealpha',0.75); hold on;
colormap(gray(100))
%%
fprintf('%s\n','use a volumetric view + slices')
FAIRfigure(3); clf; 
volView(T,omega,m,'facecolor',[240,140,100]/256,'facealpha',0.75); hold on;
colormap(gray(100))
vh = viewSlices(T,omega,m);

%%
FAIRfigure(4); clf; colordef(gcf,'black'); colormap(gray(100));
for j=1:m(3);
  vh = viewSlices(T,omega,m,'s1',64,'s2',64,'s3',[j]);
  FAIRpause(1/10);
end;
echo off
%==============================================================================
