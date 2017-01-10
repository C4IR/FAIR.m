%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
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

FAIRpause = getappdata(0,'FAIRpause');
if isempty(FAIRpause)
    FAIRpause = 'on';
end

if strcmp(FAIRpause,'on')
    builtin('pause',varargin{:});
else
%     fprintf('FAIRpause=%s\n',FAIRpause);
end
%==============================================================================
