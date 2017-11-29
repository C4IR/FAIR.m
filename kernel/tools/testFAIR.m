%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration.
% For details see
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% This is the main testing file
%==============================================================================

clear; clc;

FAIRtestStatus = getappdata(0,'FAIRtestStatus');

if isempty(FAIRtestStatus),
  FAIRtestStatus = struct(...
    'FAIRedit','off',...
    'FAIRcompile','on',...
    'FAIRrun','on',...
    'FAIRclearTemp','on'...
    );
  setappdata(0,'FAIRtestStatus',FAIRtestStatus);
  return;
end;

% clear editor files
if 0,
  allOpenFiles = matlab.desktop.editor.getAll'; % Array of open files
  fileNames = {allOpenFiles.Filename}';         % Extract file names
  K = not(strcmp([mfilename('fullpath'),'.m'],fileNames));
  allOpenFiles(K).close;
end;

% remove contents of temp folder
if strcmp(FAIRtestStatus.('FAIRclearTemp'),'on'),
  fprintf('clear folder <%s>\n',fullfile(FAIRpath,'temp'))
  cd(fullfile(FAIRpath,'temp'))
  files = dir;
  J = find(not([files(:).isdir]));
  for j=J
    delete(files(j).name);
  end;
  FAIRtestStatus.('FAIRclearTemp') = 'off'
  setappdata(0,'FAIRtestStatus',FAIRtestStatus);
  return;
end;

cd(FAIRpath);

dbstop in FAIRpause at 29

tests = {
  'testData'
  'testImgModels'
  'testTransformations'
  'testLandmarks'
  'testDistances'
  'testRegularizers'
  'testMatrixFree'
  'testNumerics'
  'testExamples'
  'testViewer'
  'testTools'
  }


for k=1:length(tests)
  edit(tests{k});
  %run(tests{k});
%   return
end;

% 
% 
% 
% 

%==============================================================================









