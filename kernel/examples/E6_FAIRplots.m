%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR:  How to use FAIRplots
%
% - loads (or generates data) based on E9_Hands_MLIR_SSD_mbElas
% - calls FAIRplotswith all options
%
% FAIRplots('clear');                  - clear the function: 
% FAIRplots('reset','mode','PIR');     - reset and set mode to parametric image registration 
% FAIRplots('init',para);              - initialize (T,R,omega,m,)
% FAIRplots('stop',para);              - display the stopping configuration
% FAIRplots('start',para);             - display the start configuration
% FAIRplots(k,para);                   - display the k'th configuration
%==============================================================================

clear; help(mfilename);

tempFile = fullfile(FAIRpath,'temp',[mfilename,'.mat'])
if ~exist(tempFile,'file'),

  % generate some dummy data to be visualized using FAIRplots
  % load data, initialize image viewer, interpolator, transformation
  E9_Hands_MLIR_SSD_mbElas
  viewImage('reset','viewImage','viewImage2D','colormap','gray(256)');
  imgModel('reset','imgModel','linearInter');
  trafo('reset','trafo','affine2D');
   
  % shortcurs for transformed image
  X = getCellCenteredGrid(omega,m);
  Rc = imgModel(dataR,omega,X);
  Tc = @(Y)  imgModel(dataT,omega,Y);
  Dc = @(Tc) SSD(Tc,Rc,omega,m);
  
  % compute the stopping, starting and final values
  builtin('clear','center')
  wStop = trafo('w0');
  Ystop = trafo(wStop,X);
  Tstop = Tc(Ystop);
  Y0    = trafo(wc,X);
  yOpt  = center(yc,m);
  T0    = Tc(Y0);
  Topt  = Tc(yOpt);
  save(tempFile,'dataT','dataR','Tc','Rc','omega','m','Tstop',...
      'Ystop','T0','Y0','Topt','yOpt')
  close all
end;
load(tempFile)

Dc = @(Tc) SSD(Tc,Rc,omega,m);
viewImage('reset','viewImage','viewImage2D','colormap',bone(256),'axis','off');
imgModel('reset','imgModel','linearInter');
trafo('reset','trafo','affine2D');

% reset FAIRplots to have a fresh start
FAIRplots('clear');

% setup the plotting functionality, no plots yet
FAIRplots('reset','mode','PIR')%mfilename);

% initialize the plots, show T and R
FAIRplots('init',struct('Tc',dataT,'Rc',dataR,'omega',omega,'m',m)); 
FAIRpause(2);

% show the stopping values
para  = struct('Tc',Tstop,'Rc',Rc,'omega',omega,'m',m,'yc',Ystop,'Jc',Dc(Tstop));
FAIRplots('stop',para);  
FAIRpause(2);

% show the starting values
para.Tc = T0; para.yc = Y0; para.Jstop = para.Jc; para.Jc = Dc(T0);
FAIRplots('start',para);   
FAIRpause(2);

% show the final values
para.Tc = Topt; para.yc = yOpt; para.Jc = Dc(Topt);
para.normdY = norm(yOpt-Y0)/norm(Ystop);
FAIRplots(1,para);
%==============================================================================
