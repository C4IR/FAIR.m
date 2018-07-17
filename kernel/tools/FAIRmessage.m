%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% Generates a nicely formatted message
%==============================================================================

function FAIRmessage(str,c)

if nargin == 0,
    FAIRmessage('nice message')
    return;
elseif (nargin == 1) && (length(str) == 1),
  fprintf('>> %s\n',char(ones(1,78))*str);
else
  if ~exist('c','var'), c = '=';  end;
  fprintf('>> %s  [ %s ]  % s\n',...
    char(ones(1,10)*c),str,char(ones(1,60-length(str))*c));
end;
%==============================================================================
