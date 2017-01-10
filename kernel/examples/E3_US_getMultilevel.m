%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: creates a multilevel data representation for an US image
%
%==============================================================================

clear, close all, help(mfilename); echo on

% load data, define a doman and an initial discretization
dataT  = double(imread('US.jpg'));
omega  = [0,size(dataT,1),0,size(dataT,2)];
m      = [128,128];

% setup image viewer
viewImage('reset','viewImage','viewImage2D','colormap','gray(256)');

% creating a multilevel representation using getMultilevel.m
ML = getMultilevel(dataT,omega,m);

% display multilevel representation of level 3
disp('ML{3}='); disp(ML{3})
echo off
%==============================================================================
