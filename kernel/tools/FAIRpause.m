%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
% ##2 
%------------------------------------------------------------------------------
% sets 'FAIRpause' on/off 
% used in the testing environment
%==============================================================================
function FAIRpause(varargin)


if nargin==1 && isstr(varargin{1}) % activate/deactivate pause status globally
    if strcmp(varargin{1},'off')
        setappdata(0,'FAIRpause','off')
        fprintf(2,'warning: pause has been disabled, see FAIRpause for details\n');        
    elseif strcmp(varargin{1},'on')
        setappdata(0,'FAIRpause','on')
        fprintf(2,'warning: pause has been enabled, see FAIRpause for details\n');        
    else
      fprintf('option <%s> is not supported, use on/off\n',varargin{1});
      error(mfilename)
    end
    return
end

FAIRtestStatus = getappdata(0,'FAIRtestStatus');

if isfield(FAIRtestStatus,'FAIRrun') ...
    & strcmp(FAIRtestStatus.('FAIRrun'),'on')
  fprintf('FAIRrun == ''on'' => no pause %s\n',mfilename);
  return;
end;

builtin('pause',varargin{:});
%==============================================================================
