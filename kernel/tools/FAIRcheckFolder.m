%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%------------------------------------------------------------------------------
% 1. locate the contents file of a folder and determine the debit
%    (list of required files), assert that contents.m is present
% 2. get a list of all files in the folder, get rid of unwanted stuff
% 3. make a balance what is missing and what is on top, list these files
% 4. sort the file list with respect to known extensions, return the list
%==============================================================================


function [folder,files,extensions] = FAIRcheckFolder(varargin)

caller  = dbstack;               % identify the name of the calling function
caller  = caller(min(length(caller),3)).name;

if nargin == 0,
  if strcmp(caller,'FAIReval')
    return;
  end;
  folder  = fileparts(which(caller));
else
  folder = varargin{1};
end;

FAIRmessage(mfilename)
fprintf('checks on [%s]\n',folder)
fprintf('caller is <%s>\n',caller)
FAIRmessage('=');

files   = {};
errmsg  = [];


% get the debits
try
  cdAct = pwd;
  cd(folder);
  debit = contents;
  cd(cdAct);
catch
  assert( 0, '''contents.m'' not executable');
end;

% get the credit
try
  credit = dir(fullfile(folder));
  credit = {credit(:).name}';
catch
  assert(0, '''dir'' not executable');
end;

% remove unwanted credits
K = [find(strcmp(credit,'.')),...
  find(strcmp(credit,'..')),...
  find(strcmp(credit,'.DS_Store')),...
  find(strcmp(credit,'.svn')),...
  find(strcmp(credit,'.git')),...
  find(strcmp(credit,'temp')),...
  find([cellfun(@(x) x(end),credit) == '~'])',...
  []];

K = unique(K);
credit(K) = [];

% check what is in or out
missing     = [];
for j=1:length(debit)
  k = find(strcmp(credit,debit{j}),1,'first');
  if isempty(k)
    missing(end+1) = j;
  else
    credit(k) = [];
  end;
end

% completeness of folder
fprintf('folder <%s> %d missing files\n',folder,length(missing));
if length(missing) > 0,
  fprintf(2,'folder <%s> %d missing files\n',folder,length(missing));
  for j=1:length(missing)
    fprintf(2,'  - %4d-of%4d missing %-30s\n',...
      j,length(missing),debit{missing(j)});
  end;
  assert( 1==0, sprintf('folder <%s> %d missing files\n',folder,length(missing)));
  return;
end;

% additional files in folder

FAIRignoreCompiles = FAIRtestPara('get','FAIRignoreCompiles')

if strcmp(FAIRignoreCompiles,'on')
  ignore = zeros(length(credit),1);
  for k=1:length(credit)
    [~,~,ext] = fileparts(credit{k});
    ignore(k) =  any(strcmp(ext,{'.o',['.',mexext]}));
  end;
  credit(find(ignore)) = [];    
end;

addons = length(credit);
% fprintf('folder <%s> %d additional files\n',folder,addons);
if addons>0,
  fprintf(2,'folder <%s> %d additional files :\n',folder,addons);
  errmsg = [errmsg,sprintf('folder <%s> %d additional files\n',folder,addons)];
  for j=1:addons
    fprintf('  - %4d-of%4d additional %-30s\n',j,addons,credit{j});
  end;
  keyboard
end;

% merge debit and andons, sort by type
files = {debit{:},credit{:}};
ext           = @(str) str(max(1,find(str=='.',1,'last')):end);
extensions    = cellfun(ext,files','UniformOutput',0)';
[filelist,OK] = sortFiles(extensions);
assert( OK, sprintf('folder <%s> contains unknown filetype(s)',folder));
files         = files(filelist);
extensions    = extensions(filelist);

FAIRmessage('=');
%------------------------------------------------------------------------------

function [C,OK] = sortFiles(list)

ext = {'.cpp','.c','.h','.o',['.',mexext],'.mat','.jpg','.m'};
C = []; R = 1:length(list);
for p=1:length(ext),
  K = find(strcmp(list,ext{p}));
  C = [C,reshape(K,1,[])];
  R = setdiff(R,K);
end;
OK = isempty(R);
for k=1:length(R)
  fprintf(' - %d-of-%d, unknown extension [%s]\n',j,length(R),list{R(k)} )
end;
%==============================================================================
