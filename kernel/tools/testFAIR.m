%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% This is the main testing file
%==============================================================================
clc;

FAIRtestPara('set',...
  'FAIRedit','off',...
  'FAIRcompile','off',...
  'FAIRignoceComiles','on',...
  'FAIRrun','on',...
  'FAIRrange','full',...
  'FAIRextension','all',...
  'FAIRfilter','',...
  'FAIRerror','on',...
  'FAIRkeyboard','off');

FAIRtestPara('disp');

%% remove contents of temp folder
if 1,
  cdAct = pwd;
  cd(fullfile(FAIRpath,'temp'))
  files = dir;
  J = find(not([files(:).isdir]));
  
  for j=J
    delete(files(j).name);
  end;
  cd(cdAct);
end;

% testData
% testImgModels
% testTransformations
% testLandmarks
% testDistances
% testRegularizers
% testMatrixFree
% testNumerics
testExamples
% testViewer
% testTools

%==============================================================================









