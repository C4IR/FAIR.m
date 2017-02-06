%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% This starts the testing environment of a folder of the toolbox
%
% Based on folder/content.m we get a list of files to be tested (also missing 
% and additional files);
% Optionally, all c-files are compiled
% Calling functions (such as testStart.m) are excluded to avoid infinite recursion
% The file list may get filter based on extension (such as "*.m") or a file range
% (such as 'full' or 1:10')
% All files are handled based on the options 
% 
% see also FAIRcheckFolder.m FAIRcheckFiles.m, FAIRmake.m, FAIRtestPara.m, testEnd.m
%
% For extented insight see 
% testFAIR.m and testTools.m
%==============================================================================

function files = testStart(folder)

caller = dbstack; % identify the name of the calling function


%% 1. get a list of files to be checked
[folder,files,extensions] = FAIRcheckFolder(folder);

%% 2. compile c files
FAIRcompile = FAIRtestPara('get','FAIRcompile');
if strcmp(FAIRcompile,'on'),
  FAIRmake(folder)
end;

%% 3. remove calling functions to avoid recusions
for j=1:length(caller),
  files(strcmp(files,[caller(j).name,'.m'])) = []; 
end;

%% 4. filter with respect to file extension
FAIRextension  = FAIRtestPara('get','FAIRextension');
if not(strcmp(FAIRextension,'all')),
  ext        = @(str) str(max(1,find(str=='.',1,'last')):end);
  extension  = cellfun(ext,files','UniformOutput',0)';
  M        = find(strcmp(extension,FAIRextension));
  files    = files(M);
end;

%% 5. filter with respect file names
FAIRfilter  = FAIRtestPara('get','FAIRfilter');
if not(isempty(FAIRfilter)) && not(strcmp(FAIRfilter,'off')),
  F     = strfind(files,FAIRfilter);
  b     = find(1-cellfun(@isempty,F));
  files = files(b);
end;

%% 5. filter with respect file range
FAIRrange = FAIRtestPara('get','FAIRrange');
if strcmp(FAIRrange,'full'),
  FAIRrange = 1:length(files);
end;

FAIRrange = intersect(FAIRrange,1:length(files));
files  = files(FAIRrange);

%% 6. execute all files, see checkFiles for details
fileOK = FAIRcheckFiles(folder,files,mfilename,'FAIRcompile','off');

J = find(fileOK == 0);
for j=J,
  fprintf(2,'-- %-2-of-%2d <%s> does not perform\n',...
    j,length(J),files{j});
end;
assert(all(fileOK == 1), 'some files are not performing');
%==============================================================================

