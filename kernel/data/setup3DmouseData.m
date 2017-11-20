%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% Cardiac gated [F-18]FDG PET data of a mouse heart acquired with the
% quadHIDAC small animal PET at the European Institute for Molecular
% Imaging (EIMI), University of Muenster, Germany.
%==============================================================================

checkDataFile
if expfileExists, return; end;

[viewer,viewPara] = viewImage('reset',...
  'viewImage','imgmontage','direction','-zyx','colormap','gray(256)');
[imgmodel,imgPara] = imgModel('reset',...
  'imgModel','splineInterMex','regularizer','moments','theta',.01)

traPara     = {'trafo','affine3D'};
disPara     = {'distance','SSD'};
regPara     = {'regularizer','mbHyperElastic','alpha',1,...
  'alphaLength',1,'alphaArea',.1,'alphaVolume',2};

load('mice3D');
dataT = double(dataT);
dataR = double(dataR);
ML = getMultilevel({dataT,dataR},omega,m,'fig',2);
% save to outfile

save(expfile,'dataT','dataR','omega','m','ML',...
  'viewPara','imgPara','traPara','disPara','regPara');

checkDataFile

xc       = getCellCenteredGrid(omega,m);
Tc       = imgModel(dataT,omega,xc);
Rc       = imgModel(dataR,omega,xc);

FAIRfigure(1,'figname',mfilename); clf;
subplot(1,2,1); viewImage(Tc,omega,m); title('template');
subplot(1,2,2); viewImage(Rc,omega,m); title('reference');
%==============================================================================
