%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function FAIRmake(file, varargin)
%
% Compile c-files or delete compiled c-files.
%
% Input:
%   file        - name of a c-file or folder
%   varargin    - variable input arguments as follows
%   clean       - remove binarys
%   cores       - compile for multi-core processos (openMP)
%   verbose     - show MEX infos
%
% Examples:
%   FAIRmake('kernel'); % compile all files of FAIR kernel
%   FAIRmake('kernel','cores',2,'verbose',true); % compile all for 2 cores
%   FAIRmake('imgModels','clean',true); % delete compiled interpolation files
%   FAIRmake('linearInterMexC.cpp','clean',true); % delete compiled file
%   FAIRmake('linearInterMexC.cpp'); % compile 'linearInterMexC.cpp'
%==============================================================================

function FAIRmake(file,varargin)

if nargin == 0,
  help(mfilename)
  FAIRmake('data')
  return;
end;

overwrite = 0;     % if overwrite == 0, do not re-compile
cores     = 1;     % choose number of cores; if cores>1 openMP is used
verbose   = false; % some output

for k=1:2:length(varargin) % overwrites default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end

if verbose > 0,
    fprintf('%s aims to compile <%s>\n',mfilename,file)
end;

folder = {
    'data'
    'imgModels'
    'transformations'
    'landmarks'
    'distances'
    'regularizers'
    };
        
Make = @(file) make(file, overwrite, cores, verbose);
if any(strcmp(file,folder))
    file = fullfile(FAIRpath,'kernel',file);
end

if exist(file) == 7, % is a folder
    shortname = file(max([0,find(file==filesep,1,'last')])+1:end);

    if strcmp(shortname,'kernel') % build all files in kernel
        for k=1:length(folder)
            longname = fullfile(FAIRpath,'kernel',folder{k});
            fprintf('make C files of %s\n',longname);
            FAIRmake(longname,varargin{:})
        end;
        return
   end
    
    files = dir(file);
    files = {files(:).name}';
    ext        = @(str) str(max(1,find(str=='.',1,'last')):end);
    extensions = cellfun(ext,files','UniformOutput',0)';
    J = find( strcmp(extensions,'.cpp') | strcmp(extensions,'.c') );
    fprintf('[%s] built %d files:\n',mfilename,length(J));
    for j=1:length(J),
      Make(files{J(j)});
    end;
elseif exist(file,'file')
    Make(file);
end

%------------------------------------------------------------------------------

function make( filename, overwrite, cores, verbose )

fprintf('mex %s\n',filename)

[folder,file,ext] = fileparts(which(filename));
mexfilename = fullfile(folder,[file,'.',mexext]);

if (exist(mexfilename) == 3) & ~overwrite,
  fprintf(' -- file <%s> already exists\n',mexfilename);
  return
end;

[~, maxsize] = computer;    % determine if 64bit (maxsize==2^48-1) or 32bit
if maxsize==2^48-1,
  options = '-largeArrayDims';
else
  options = '';
end
if verbose>0,
  options = [options ' -silent'];
end

cdAct = pwd;
cd(folder)

if strcmp(filename,'geometryC.cpp'),
  cd(cdAct);
  return;
elseif strcmp(filename,'geometryMexC.cpp'),
  filename = [filename,' geometryC.cpp'];
end

if cores == 1,
  str = ['mex ' filename ' -g  CFLAGS="\$CFLAGS -p ',...
    '-ftree-vectorize " CXXFLAGS="\$CXXFLAGS -p ',...
    '-ftree-vectorize "  LDFLAGS="\$LDFLAGS -p ',...
    '-ftree-vectorize " ' options];
else
  setenv('OMP_NUM_THREADS',num2str(cores))
  str = ['mex ' filename ' -O CC=gcc CXX=g++ LD=g++ CFLAGS="\$CFLAGS ',...
    '-fopenmp -ftree-vectorize" CXXFLAGS="\$CXXFLAGS -fopenmp ',...
    '-ftree-vectorize" LDFLAGS="\$LDFLAGS ',...
    '-fopenmp -ftree-vectorize" COMPFLAGS="$COMPFLAGS -fopenmp" ',...
    options];
end;

if verbose == 1,
  fprintf('    mex %s\n',filename)
elseif verbose>1, 
  disp(str); 
end;

try
  eval(str);
catch err,
  warning(sprintf('build of mex-file <%s> was NOT successful!',filename));
end;
cd(cdAct)
%==============================================================================
