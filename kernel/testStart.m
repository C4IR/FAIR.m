%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
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

function testStart


caller        = getappdata(0,'caller');
folder        = fileparts(which(caller));
FAIRcheckList = getappdata(0,'FAIRcheckList');

if isempty(FAIRcheckList),
     FAIRcheckList = FAIRcheckFolder(folder);
     shortname = @(x) x(max(0,find(x == filesep,1,'last'))+1:end-2);
     files = cellfun(shortname,{FAIRcheckList(:).name},'UniformOutput',0);
     C     = find(strcmp(files,caller));
     FAIRcheckList(C).check = 2;
     setappdata(0,'FAIRcheckList',FAIRcheckList);
end;

FAIRcompile = FAIRtestPara('get','FAIRcompile');
if strcmp(FAIRcompile,'on'),
  FAIRmake(folder)
end;

FAIRcheckFiles

% 
% fprintf('check %d files\n',length(files))
% for j=1:length(files),
%   fprintf('%3d-of%d: %s\n',j,length(files),files{j});
% end;
%   
%   
% %% 2. compile c files
% 
% %% 3. remove calling functions to avoid recusions
% files(strcmp(files,[caller,'.m'])) = [];
% 
% %% 4. filter with respect to file extension
% FAIRextension  = FAIRtestPara('get','FAIRextension');
% if not(strcmp(FAIRextension,'all')),
%   ext        = @(str) str(max(1,find(str=='.',1,'last')):end);
%   extension  = cellfun(ext,files','UniformOutput',0)';
%   M        = find(strcmp(extension,FAIRextension));
%   files    = files(M);
% end;
% 
% %% 5. filter with respect file names
% FAIRfilter  = FAIRtestPara('get','FAIRfilter');
% if not(isempty(FAIRfilter)) && not(strcmp(FAIRfilter,'off')),
%   F     = strfind(files,FAIRfilter);
%   b     = find(1-cellfun(@isempty,F));
%   files = files(b);
% end;
% 
% %% 5. filter with respect file range
% FAIRrange = FAIRtestPara('get','FAIRrange');
% if strcmp(FAIRrange,'full'),
%   FAIRrange = 1:length(files);
% end;
% 
% FAIRrange = intersect(FAIRrange,1:length(files));
% files  = files(FAIRrange);
% 
% %% 6. execute all files, see checkFiles for details
% fileOK = FAIRcheckFiles(folder,files,mfilename,'FAIRcompile','off');
% 
%==============================================================================

