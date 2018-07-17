%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% MRI data of a human knee, image courtesy by
% Thomas Netsch, Philips Medical Solutions, Hamburg
%==============================================================================

checkDataFile
if expfileExists, return; end;

infile = fullfile(FAIRpath,'kernel','data','knee3D.mat');
try
  load(infile);
catch
  fprintf('sorry, this data is not available');
  expfile = [];
  return;
end;

dataT = double(dataT);
dataR = double(dataR);
viewK = @(T) volView(T,omega,m,'isovalue',15,'view',[-20,-5],'colormap','bone(256)');

FAIRfigure(1,'color','w','figname',sprintf('%s/template',mfilename)); clf;
viewK(dataT); hold on; axis off; colormap(gray(128));
FAIRfigure(2,'color','w','figname',sprintf('%s/reference',mfilename)); clf;
viewK(dataR); hold on; axis off; colormap(gray(128));

[viewer,viewPara] = viewImage('reset','viewImage','imgmontage',...
  'direction','-zyx','colormap','bone(256)');
imgPara     = {'inter','linearInterMex'};
traPara     = {'trafo','affine3Dsparse'};
disPara     = {'distance','SSD'};
regPara     = {'regularizer','mfElastic','alpha',500,'mu',1,'lambda',0};

ML = getMultilevel({dataT,dataR},omega,m,'fig',2);

% save to outfile
save(expfile,'dataT','dataR','omega','m','ML',...
  'viewPara','imgPara','traPara','disPara','regPara');

checkDataFile
%==============================================================================
