%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% This files deals the parameters for the Testing environment
%
% options:
% FAIRtestPara('disp')          displays the current settings
% FAIRtestPara('set',list{:})	sets defaults and overwrites with user supplied values
% FAIRtestPara('get',name)		gets value of variable name
% FAIRtestPara('clear')			clears all variables
%
%==============================================================================

function varargout = FAIRtestPara(varargin)

caller = dbstack; % identify the name of the calling function

if nargin == 0,
  help(mfilename);
  FAIRtestPara('disp');
  return;
end;

if nargin > 0 && isstr(varargin{1}),
  task = varargin{1};
  varargin(1) = [];
  
  switch task,
    case 'set',
      setFAIRtestpara(varargin{:});

    case 'get',     
      for j = 1:length(varargin)
        varargout{j} = getappdata(0,varargin{2*j-1});
      end;

%    case 'value',
%      varargout{1} = strcmp(getappdata(0,varargin{1}),'on');
 
    case 'disp',
      dispFAIRtestpara;
      
    case 'clear',
      FAIRpara = getappdata(0);
      Fnames   = fieldnames(FAIRpara);
      K        = reshape(find(cellfun(@isempty,strfind(Fnames,'FAIR'))==0),1,[]);
      for k = K
        rmappdata(0,Fnames{k});
      end;
    otherwise, error('1');
  end;
else
  error(1)
end;

%------------------------------------------------------------------------------

function setFAIRtestpara(varargin)

para = {
  'FAIRfolder',         ''
  'FAIRedit',           'on'
  'FAIRcompile',        'on'
  'FAIRignoreCompiles', 'on'
  'FAIRrun',            'on'
  'FAIRrange',          'full'
  'FAIRextension',      'all'
  'FAIRerror',          'on'
  'FAIRkeyboard',       'on'
  'FAIRfilter',         'off'
  'FAIRpause',          'off'
  'FAIRinput',          'off'
};

for j=1:size(para,1)
  value = getappdata(0,para{j,1});
  if isempty(value)
    setappdata(0,para{j,1},para{j,2})
  end;
end;

for k=1:length(varargin)/2, % overwrites default parameter
    setappdata(0,varargin{2*k-1},varargin{2*k});
end;

%------------------------------------------------------------------------------

function dispFAIRtestpara;
FAIRpara = getappdata(0);
Fnames   = fieldnames(FAIRpara);
K        = reshape(find(cellfun(@isempty,strfind(Fnames,'FAIR'))==0),1,[]);

if isempty(K)
return;
end;

FAIRmessage('.');
for k = K
name  = Fnames{k};
value = getappdata(0,name);
if ischar(value)
  fprintf('%-20s: %s\n',name,value)
elseif isnumeric(value) && (length(value) == 1)
  fprintf('%-20s= %s\n',name,num2str(value))
elseif isnumeric(value),
  fprintf('%-20s: %d files\n',name,length(value))
else
  keyboard
end;
end;
FAIRmessage('.');

%==============================================================================
