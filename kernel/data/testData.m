%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% This is a testing environment for the files in the folder kernel/data
% 1. Based on data/contents, a list of required files is generated and it 
%    is verified, that all files are present; additional files are listed.
% 2. All c-files are compiled.
% 3. All files are executed.
% 4. Check if all expected variables are present.
%==============================================================================


FAIRcheckFiles(mfilename);

FAIRtestStatus = getappdata(0,'FAIRtestStatus')

shortname = @(x) x(max(0,find(x == filesep,1,'last'))+1:end);
files = cellfun(shortname,{FAIRtestStatus.('testData')(:).name},'UniformOutput',0);

%% Check if all expected variables are present
folder        = fileparts(which(mfilename));
J = find(cellfun(@(s) isequal(strfind(s,'setup'),1), files));
fprintf('folder <%s> checking %d .mat-files\n',folder,length(J));

mfiles = cellfun( @(s) s(1:end-2), files(J),'UniformOutput',0);
credit   = {'dataT','dataR','omega','omega','m','ML','viewPara'};

for j = 1:length(mfiles)
  fprintf(2,'  %4d-of%4d check contents of %-30s\n',j,length(J),mfiles{j});
  
  % check file contents
  debit = whos('-file',mfiles{j});
  debit = {debit(:).name};
  OK = 1;
  for k=1:length(credit),
    if isempty(find(strcmp(debit,credit{k})==1));
      fprintf(2,' variable [%s] is missing!\n',credit{k})
      OK = 0;
    end;
  end;
  assert( OK == 1, sprintf('<%s> data is incomplete',mfiles{j}))
end;

testEnd;
%==============================================================================
