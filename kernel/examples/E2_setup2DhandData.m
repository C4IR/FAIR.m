%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
% ##2
%==============================================================================
% Template for setup data
%     (here: use hand xray's; see setup2DhandsData.m for references)
%
% This file initializes the following data (which can then be saved): 
%   dataT      template  image, a d-array of size mD, 
%   dataR      reference image, a d-array of size nD, 
%   omega      domain specification
%              omega = (omega(1),omega(2)) x  (omega(3),omega(24))
%   m          initial discretization 
%   ML         multi-level representation of the data
%   L M        landmarks, if available
%
% For representation and visualization
%   viewPara   options for image viewer
%   imgPara    options for image interpolation
% see also setup2DhandsData and E2_setupBrainData
%==============================================================================

clear, close all, help(mfilename), 
echo on


% note: FAIR uses a right handed x-y coordinate system
%   in order to convert the image i-j coordinated we use 
%   flipud()  and transpose() converts ij to xy coordinates
%   double    converts from uint to double
%   rgb2gray  converts from rgb to gray scale

Tdata = double(flipud(imread('hands-T.jpg'))');
Rdata = double(flipud(imread('hands-R.jpg'))');
mD    = size(Tdata);
nD    = size(Tdata);

omega = [0,20,0,25];    % specify physical domain
m     = [128,256];      % note the resolution of the representation might be 
                        % different to the data size     

% for this data landmarks (LM) have been identified
% note: LM are physical (x-y coordinates) in the omega domain

LM = [
   5.5841   17.2664    2.6807   12.7797
  10.7243   21.6121    7.2028   19.6795
  13.2477   21.6121   10.1865   20.8916
  15.2570   19.2290   12.5175   20.0991
  15.8645   15.1636   14.3357   16.7424
   5.3972    8.1075    7.9953    6.3462
   7.5000    5.9579   11.8648    5.6469
  ];

% set view options and interpolation options and initialize viewer and interpolator
viewPara = {'viewImage','viewImage2D','colormap','gray(256)'};
viewImage('reset',viewPara{:});

imgPara = {'imgModel','linearInter'};
imgModel('reset',imgPara{:});

% create multilevel representation of the data
ML = getMultilevel({Tdata,Rdata},omega,m,'fig',2);

% visualize
xc = getCellCenteredGrid(omega,m);

FAIRfigure(1,'figname',mfilename); clf;
  subplot(1,2,1); viewImage(imgModel(Tdata,omega,xc),omega,m); hold on;
    ph = plotLM(LM(:,1:2),'numbering','on','color','r'); 
    set(ph,'linewidth',2,'markersize',20);
    title(sprintf('%s','template&LM'),'fontsize',20);
  subplot(1,2,2); viewImage(imgModel(Rdata,omega,xc),omega,m); hold on;
    ph = plotLM(LM(:,3:4),'numbering','on','color','g','marker','+');
    set(ph,'linewidth',2,'markersize',20);
    title(sprintf('%s','reference&LM'),'fontsize',20);
    
echo off
%==============================================================================
