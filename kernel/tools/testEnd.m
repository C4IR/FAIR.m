%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% finishes a test run, nice output, small break, close all
%==============================================================================

function testEnd

caller = dbstack; % identify the name of the calling function
caller = caller(min(length(caller),2)).name;
caller = which(caller);

str = {
  sprintf('%% %s\n',char(ones(1,78)*'#'))
  sprintf('    file <%s> terminated\n',caller)
  sprintf('%% %s\n',char(ones(1,78)*'#'))
};

fprintf('%s\n',str{:});

builtin('pause',2);
close all

%==============================================================================
