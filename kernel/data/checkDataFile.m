%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% This scripts checks whether the required data is already available.
% If if so and possible, it also initializes the modules
% viewimage, imgModel, trafo, distance, regularizer for later use;
% see setup2DhandsData.m for an example
%==============================================================================
caller          = dbstack;        % identify the name of the calling function
caller          = caller(min(length(caller),2)).name;
expfile         = fullfile(FAIRpath,'temp',[caller,'.mat']);
expfileExists   = (exist(expfile,'file') == 2);

if ~expfileExists,
  return;         % matfile has to be generated
end;

fprintf('%s: load(%s); and initialize modules\n',...
  caller,expfile(length(FAIRpath)+1:end));
clear caller
load(expfile);

% initialize toolbox
if ~isempty(who(expfile,'viewPara')),    viewImage('reset',viewPara{:});  end;
if ~isempty(who(expfile,'imgPara')),     imgModel('reset',imgPara{:});    end;
if ~isempty(who(expfile,'traPara')),     trafo('reset',traPara{:});       end;
if ~isempty(who(expfile,'disPara')),     distance('reset',disPara{:});    end;
if ~isempty(who(expfile,'regPara')),     regularizer('reset',regPara{:}); end;
%==============================================================================
