%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: density estimation using a histogram
%
%==============================================================================

clear, close all, help(mfilename);

n = 1000; 
x = linspace(0,2*pi,n); T = 0.5*(sin(x)+1); % discretized function
minT = 0; maxT = 1;                         % bounds for the bins
numberBins = 5;                             % number of bins
binWidth   = (maxT-minT)/numberBins;        % bin width
bins       = 0:binWidth:maxT;               % the bins
binsExt    = [-inf,bins(2:end-1),inf];      % don't miss anything

rhoHat = histc(T,binsExt); % compute the histogram and plot it
bar(bins+binWidth/2,rhoHat,0.99,'edgecolor','w','facecolor',0.8*[1,1,1]);
axis([minT,maxT,0,inf])

%==============================================================================
