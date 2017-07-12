% ==================================================================================
% (c) Andreas Mang 2016/09/23, see FAIR.2 and FAIRcopyright.m.
% (http://www.mic.uni-luebeck.de/)
dataset='na06';
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
max(dataT(:))

%omega = [0,m(1),0,m(2),0,m(3)];
omega = [0,20,0,23.4375,0,20];
ML = getMultilevel({dataT,dataR},omega,m,'fig',2);
save(outfile,'dataT','dataR','dataTl','dataRl','omega','m','MLdata');
save(outfile,'-append','viewOptn','imgOptn','traOptn','disOptn','regOptn');
checkSetupDataFile;

