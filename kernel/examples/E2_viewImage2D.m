% ==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% 
% Tutorial for FAIR: shows how to visulaize 2D data with FAIR.
%
% - load data             (here 'US.jpg')
% - setup viewer          (viewImage2D)
% - view  data
%==============================================================================

clear, close all, help(mfilename);
echo on

% load data, get size and specify physical domain
dataT = double(imread('US.jpg'));
m     = size(dataT);
omega = [0,m(1),0,m(2)]; 

% setup image viewer
viewImage('reset','viewImage','viewImage2D','colormap','gray(256)','axis','off');

% view data
viewImage(dataT,omega,m,'title','my first image');

echo off
%==============================================================================
