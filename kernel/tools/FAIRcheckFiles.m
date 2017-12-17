%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration.
% For details see
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration.
% For details see
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% Evaluates a list of files respecting the file type
% parameterized via FAIRtestPara.m
%
%
%==============================================================================

function FAIRcheckFiles(caller)

if nargin == 0,
  FAIRtestStatus = getappdata(0,'FAIRtestStatus');
  caller = FAIRtestStatus.('caller');
end;

FAIRtestStatus = getappdata(0,'FAIRtestStatus');
FAIRtestStatus.('caller') = caller;
FAIRtestStatus.('folder') = fileparts(which(caller))

if not(isfield(FAIRtestStatus,caller))
  FAIRtestStatus.(caller) = {};
end;

if isempty(FAIRtestStatus.(caller)),
  FAIRtestStatus.(caller) = FAIRcheckFolder(FAIRtestStatus.('folder'));
  shortname = @(x) x(max(0,find(x == filesep,1,'last'))+1:end-2);
  filename  = cellfun(shortname,{FAIRtestStatus.(caller)(:).name},'UniformOutput',0);
  
  excludes = {...
    caller,...
    'testFAIR',...
    'FAIRcheckFolder',...
    'FAIRcheckFiles',...
    'testEnd',...
    'getLandmarks'};
  
  for j=1:length(excludes)
    C     = find(strcmp(filename,excludes{j}));
    if ~isempty(C),
      FAIRtestStatus.(caller)(C).check = 2;
    end;
  end,
  setappdata(0,'FAIRtestStatus',FAIRtestStatus);
end;

if strcmp(FAIRtestStatus.('FAIRcompile'),'on'),
  FAIRmake(FAIRtestStatus.('folder'))
end;

%------------------------------------------------------------------------------
% main loop, run over all files
%------------------------------------------------------------------------------
mexall = mexext('all');
for j=1:length(FAIRtestStatus.(caller))
  name  = FAIRtestStatus.(caller)(j).name;
  check = FAIRtestStatus.(caller)(j).check;
  [~,file,ext] = fileparts(name);
  
  fprintf('%3d-of-%d %-40s : check = %3d\n',...
    j,length(FAIRtestStatus.(caller)),[file,ext],check);
  
  if check == 1 | check == 2 | check == 3,
    builtin('pause',0.1);
  else
    file = [file,ext];
    
    switch ext,
      case '.mat',    OK = FAIRload(file);
      case '.jpg',    OK = FAIRjpg(file);
      case '.m',
        OK = 1;
        if strcmp(FAIRtestStatus.('FAIRrun'),'on'),
          OK = FAIReval(name);
        end;
        if strcmp(FAIRtestStatus.('FAIRedit'),'on') || OK==0,
          FAIRopen(name);
        end;
      case {'.cpp','.c'},
        if strcmp(FAIRtestStatus.('FAIRcompile'),'on'),
          OK = FAIRbuild(file);
        else
          OK = 1;
        end;
      case {'.h','.o'},
        OK = 1;

      case strcat({'.'}, {mexall.ext}),
        OK = 3;

        
      otherwise,
        ext
        error('12');
    end;
    
    FAIRtestStatus.(caller)(j).check = OK;
    setappdata(0,'FAIRtestStatus',FAIRtestStatus);
  end;
end;

J = find([FAIRtestStatus.(caller)(:).check] <= 0);

if length(J)>0,
  fprintf('non-performing files:\n');
  for j=J
    fprintf('  - %4d-of-%4d %-30s\n',...
      j,length(FAIRtestStatus.(caller)),FAIRtestStatus.(caller)(j).name);
  end;
end;

%------------------------------------------------------------------------------
function OK = FAIRload(file)

try
  load(file)
  OK = 1;
catch
  OK = -1;
end;
%------------------------------------------------------------------------------
function OK = FAIRjpg(file)
try
  B = imread(file);
  figure(1); clf;
  set(1,'position',[3400 900 500 500]);
  imagesc(B); axis image
  builtin('pause',1)
  OK = 1;
catch
  OK = -1;
end;
%------------------------------------------------------------------------------

function OK = FAIReval(file)
% clear = @()         fprintf(2,'clear has been disabled for auto-testing\n');
file = file(1:end-2);
try
  eval(sprintf('dbclear in %s',file))
  run(file);
  OK = 1;
catch e
  fprintf('file does not run without errors\n');
  OK = -10;
end;
%------------------------------------------------------------------------------
function OK = FAIRbuild(file)
try
  FAIRmake(file)
  OK = 1;
catch
  OK = -1;
end;
%------------------------------------------------------------------------------
function OK = FAIRopen(file)
try
  edit(file)
  OK = 1;
catch
  OK = -1;
end
%==============================================================================
