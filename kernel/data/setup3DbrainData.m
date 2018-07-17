%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% MRI data of a human brain,
% data courtesy Ron Kikinis, Brigham & Women's Hospital, Boston, USA
%==============================================================================

checkDataFile
if expfileExists, return; end;

load('brain3D')

% setup image viewer,
[viewer,viewPara] = viewImage('reset','viewImage','imgmontage',...
  'colormap','gray(256)','direction','-zyx');

FAIRfigure(1); clf;
subplot(1,2,1); viewImage(dataT,omega,m);
subplot(1,2,2); viewImage(dataR,omega,m);

ML = getMultilevel({dataT,dataR},omega,m,'fig',2);

imgPara  = {'imgModel','splineInterMex'};
traPara  = {'trafo','affine3Dsparse'};
disPara  = {'distance','SSD'};
regPara  = {'regularizer','mfElastic','alpha',1e3,'mu',1,'lambda',0};

% save to outfile
save(expfile,'dataT','dataR','omega','m',...
  'ML','viewPara','imgPara','traPara','disPara','regPara');
checkDataFile
%==============================================================================
