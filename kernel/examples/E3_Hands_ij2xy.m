%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: demo for 2D data setup and visualization
%
% load images                        (hands-?.jpg)
% visualize ij and omega based data  (viewImage)
% creates multi-lebel representation (getMultilevel) 
% 
%==============================================================================

clear, close all, help(mfilename); 

% load data
Tij = imread('hands-T.jpg'); Tdata = double(flipud(Tij))';
Rij = imread('hands-R.jpg'); Rdata = double(flipud(Rij))';

% specify domain Omega = [omega(1),omega(2)]x[omega(3),omega(3)];
omega = [0,20,0,25];
m     = size(Tdata);

% setup image viewer
viewImage('reset','viewImage','viewImage2D','colormap','gray(256)');

% visualize
FAIRfigure(1);  colormap(gray(256));
subplot(2,2,1); imagesc(Tij);              title('original T data, uint8, ij');
subplot(2,2,2); viewImage(Tdata,omega,m);  title('FAIR T data on \Omega, xy');
subplot(2,2,3); imagesc(Rij);              title('original R data, uint8, ij');
subplot(2,2,4); viewImage(Rdata,omega,m);  title('FAIR R data on \Omega, xy');

% create multi-level representaion
ML = getMultilevel({Tdata,Rdata},omega,m,'fig',2);

%==============================================================================
