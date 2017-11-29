%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
% ##2 
%------------------------------------------------------------------------------
% sets 'FAIRinput' on/off 
% used in the testing environment
%==============================================================================

function status = FAIRinput(varargin)

status = getappdata(0,'FAIRinput');
if isempty(status),
  status = 'on';
  setappdata(0,'FAIRinput',status);
end;

if nargin == 0,
  return;
elseif nargin == 1,
  status = varargin{1};
  setappdata(0,'FAIRinput',status);
  return;
end;

fprintf('%%----- %-24s %-30s',[varargin{1},':'],varargin{2});
if strcmp(status,'on'),
  status = input('<0,1> : ');
else
  status = varargin{3}
end;

%==============================================================================
    
    