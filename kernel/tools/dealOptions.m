%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% This functions deals persistent variables in administrative functions such as
% viewImage, viewImage, imgModel, trafo, distance, regularizer, ....
%
% Input:
%   optn        a list of current options, to be updated
%               cell or struct, 2*n entries, e.g.
%               optn.omega = [1,2];
%               optn.m     = 14;     
%               or optn    = {'omega',[1,2],'m',14;}
%   varargin    a list of new options, cell with 2*p entries
%
% Output:
%   method      specifies the particular method to be used in caller, see example
%   optn        updated list of options
%   task        task to be performed either in options or in caller
%   stop        flag for caller to terminate 
%
% Example:
% caller='imgModel', method='splineInter', p1='theta', pV1=10
% 
% The followings modes can be uses
%    1.    nargin == 0,
%         help, run minimal example, return
%

% elseif    ~ischar(varargin{1}) -> task='none'; stop = 0;
% else      task = varargin{1} end one of the following
%
% task = 'reset', caller('reset',caller,method,p1,pV1,...), stop = 1
%   set optn={caller,method,p1,pV1,...}
%
% task = 'set',   caller('set',x1,xV1,...),                 stop = 1
%   update optn by overwiting (if x1=pk) or adding to list
%
% task = 'get',   caller('get',optn,pk),                    stop = 1
%   if field exists, method = OPTN.(pk), else method = [];
%
% task = 'disp',  caller('disp'),                           stop = 1
%   displays all the optn, fieldname:optn.(fieldname)
%
% task = 'clear', caller('clear'),                          stop = 1
%   set optn = [];
%
%  otherwise: -> task, stop = 0

%==============================================================================


function [method,optn,task,stop] = dealOptions(caller,optn,varargin)

method = [];                    % initialize method
% caller = dbstack;               % identify the name of the calling function
% caller = caller(min(length(caller),2)).name;

task   = [];
stop   = 1;

if nargin == 0,  
  help(mfilename);
  runMinimalExample;
  return;
end;

% transform optn from cell to struct
if iscell(optn),
  optn = cell2struct(optn(2:2:end),optn(1:2:end),2);
end;

% check whether work is to be done here, i.e. the first input argument 
% varargin{1} is a str
if nargin>2 && ischar(varargin{1})

  task        = varargin{1};
  varargin(1) = [];

  switch task
    case {'reset','set'},       % reset/set options
      if strcmp(task,'reset'),  % clear optn
        optn = []; 
      end;   

      % update the field with the value from the varargin list
      for k=1:2:length(varargin),
        optn.(varargin{k}) = varargin{k+1};
      end;

    case 'disp',                % display the options
      fprintf('<%s> options\n',caller); 
      if isempty(optn), 
        fprintf('is empty\n'); 
      else
        disp(optn)
      end;

    case 'clear';               % clear the options
      optn   = []; 
      fprintf('--- cleared OPTN in <%s> ---\n',caller); 

    case 'get',
      for k=1:length(varargin),
        method = getOption(optn,varargin{1});
      end;
      return
    
    otherwise, 
      stop = 0;
    
  end;

  [method,optn] = getOption(optn,caller);
  return;

else

  [method,optn] = getOption(optn,caller);
  stop = (nargin < 3);
  return;

end;


%----------------------------------------------------------------------------------------

function [value,optn] = getOption(optn,field);

% initialize value= []; update by the value in optn
value = []; if isfield(optn,field),  value = optn.(field); end;
% transfer struct optn to list
if isstruct(optn),
  fields = fieldnames(optn);
  values = struct2cell(optn);
  optn = cell(1,2*length(fields));
  optn(1:2:2*length(fields)-1) = fields;
  optn(2:2:2*length(fields))   = values;
end;

%------------------------------------------------------------------------------

function runMinimalExample
% initialize  parameters
imgModel('reset','imgModel','splineInter','regularizer','moments','theta',1);

% show the configuration
imgModel('disp');

% get scheme and its parameters
[scheme,parameter] = imgModel;

% clear the parameters
imgModel('clear');
%==============================================================================

