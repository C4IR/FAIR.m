% =========================================================================
% This code is part of the Matlab-based toolbox 
% LagLDDDM - A Lagrangian Gauss--Newton--Krylov Solver for Mass- and 
%                        Intensity-Preserving Diffeomorphic Image Registration
% 
% For details and license info see 
% - https://github.com/C4IR/FAIR.m/tree/master/add-ons/LagLDDMM
%
% Setup NIREP problem (data needs to be obtained 
% separately) 
%
% =========================================================================
dataset='na02';
example = ['3D-nirep-',dataset];
checkSetupDataFile; if OK, return; end;

% set view options and interpolation options
[viewer,viewOptn] = viewImage('reset','viewImage','imgmontage','colormap','gray(256)','direction','-zyx');

% setup interpolation scheme
imgOptn  = {'imgModel','splineInterMex','regularizer','moments','theta',.01};

% setup transformation used in the parametric part
traOptn  = {'trafo','affine3D'};

% setup distance measure
disOptn  = {'distance','SSD'};

% initialize the regularizer for the non-parametric part
regOptn = {'regularizer','mbHyperElastic','alpha',1,'alphaLength',1,'alphaArea',.1,'alphaVolume',2};

FAIRmessage(mfilename)
load('na01-128x150x128');
load('na01-128x150x128-labels');
load([dataset,'-128x150x128']);
load([dataset,'-128x150x128-labels']);
dataT = double(dataT);
dataR = double(dataR);
dataTl = double(dataTl);
dataRl = double(dataRl);
m = size(dataT);

dataT = 256.*dataT;
dataR = 256.*dataR;

%omega = [0,m(1),0,m(2),0,m(3)];
omega = [0,20,0,23.4375,0,20];
ML = getMultilevel({dataT,dataR},omega,m,'fig',0);
%ML = getMultilevel({dataT,dataR},omega,m);
%save(outfile,'dataT','dataR','omega','m','ML');
save(outfile,'dataT','dataR','dataTl','dataRl','omega','m','ML');
save(outfile,'-append','viewOptn','imgOptn','traOptn','disOptn','regOptn');
checkSetupDataFile;

% xc       = getCellCenteredGrid(omega,m);
% viewData = @(I) viewImage(imgModel(I,omega,xc),omega,m);

% FAIRfigure(1,'figname',mfilename); clf;
% subplot(1,2,1); viewData(dataT); title('template');
% subplot(1,2,2); viewData(dataR); title('reference');
