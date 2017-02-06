%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% Evaluates a list of files respecting the file type
% parameterized via FAIRtestPara.m
%
%
%==============================================================================

function fileOK =checkFiles(folder,files,caller,varargin)

if nargin == 0,
  return;
end;

FAIRtestPara('set',varargin{:})
FAIRtestPara('disp');
FAIRedit     = FAIRtestPara('get','FAIRedit');
FAIRcompile  = FAIRtestPara('get','FAIRcompile');
FAIRrun      = FAIRtestPara('get','FAIRrun');
FAIRkeyboard = FAIRtestPara('get','FAIRkeyboard');
FAIRerror    = FAIRtestPara('get','FAIRerror');

fprintf('%s: folder = <%s>\n',mfilename,folder);
fprintf('%-20s: %s\n','range',sprintf('%d files',length(files)));

% dusplay files to be processed
for j=1:length(files),
  fprintf('  - %4d-of-%4d %-30s\n',j,length(files),files{j});
end;

%------------------------------------------------------------------------------
% main loop, run over all files
%------------------------------------------------------------------------------
fileOK = ones(length(files),1);
for j=1:length(files),

  file = files{j};
  ext  = file(max(1,find(file=='.',1,'last')):end);
  
  FAIRmessage('=');
  fprintf(2,'  - %4d-of-%4d %-30s\n',j,length(files),files{j});
  FAIRmessage('=');
  
  if strcmp([caller,'.m'],file),
    ext = 1;
  end;
  
  switch ext,
    case 1,
      fprintf('do not execute the caller <%s>\n',file)
      OK = 1;
    case '.mat',    OK = FAIRload(file);
    case '.jpg',    OK = FAIRjpg(file);
    case '.m',  
      OK = 1;
      if strcmp(FAIRrun,'on'),
        OK = FAIReval(fullfile(folder,file));
      end;
      if strcmp(FAIRedit,'on') || not(OK),
        FAIRopen(fullfile(folder,file));
      end;
    case {'.cpp','.c'},
      if strcmp(FAIRcompile,'on'),
        OK = FAIRbuild(file);
      else
        OK = 1;
      end;
    case {'.h','.o',lower(['.',mexext])},
      OK = 1;
      
    otherwise,
      ext
      keyboard
      error('unkown filetype');
  end;
  
  fileOK(j) = OK;
  if not(OK)
    fprintf(2,'file <%s> does not perform and error mode is [%s]\n',...
      file,FAIRerror);
    if strcmp(FAIRerror,'on'),
      error('break')
    end;
  end;
%   assert( OK == 1, ...
%     sprintf('file  %4d-of-%4d <%s> has problems',j,length(files),file) );
  if ~OK,
    sprintf('file  %4d-of-%4d <%s> has problems',j,length(files),file);
  end;
  
  if strcmp(FAIRkeyboard,'on'),
    fprintf('[keyboard=''on'', see testStart for options]\n')
    keyboard;
  end;
  
end;
%------------------------------------------------------------------------------


% display performing files
J = find(fileOK == 1);
if length(J) > 0,
  fprintf(2,'performing files:\n')
  for j=1:length(J)
    fprintf('  - %4d-of-%4d %-30s\n',j,length(J),files{J(j)});
  end;
end;

% display non-performing files
J = find(fileOK == 0);
if length(J) > 0,
  fprintf(2,'non-performing files:\n')
  for j=1:length(J)
    fprintf(2,'  - %4d-of-%4d %-30s\n',j,length(J),files{J(j)});
  end;
end;
FAIRmessage('=');
builtin('pause',1)

%------------------------------------------------------------------------------
function OK = FAIRload(file)

try
  load(file)
  OK = 1;
catch
  OK = 0;
end;
%------------------------------------------------------------------------------
function OK = FAIRjpg(file)
try
  B = imread(file);
  figure(1); 
  set(1,'position',[3400 900 500 500]);
clf;
  imagesc(B); axis image
  builtin('pause',1)
  OK = 1;
catch
  OK = 0;
end;
%------------------------------------------------------------------------------

function OK = FAIReval(file)
% clear = @()         fprintf(2,'clear has been disabled for auto-testing\n');
file = file(1:end-2);
try
  run(file);
  OK = 1;
catch
  fprintf('file does not run without errors\n');
  OK = 0;
end;
%------------------------------------------------------------------------------
function OK = FAIRbuild(file)
try
  FAIRmake(file)
  OK = 1;
catch
  OK = 0;
end;
%------------------------------------------------------------------------------
function OK = FAIRopen(file)
try
  edit(file)
  OK = 1;
catch
  OK = 0;
end
%==============================================================================
