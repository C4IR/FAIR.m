%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function D = sdiag(d)
%
% returns sparse diagonal matrix 
%
%     | d(1)              |
%     |     d(2)          |
% D = |          \        | in R^{n,n},   n=numel(d)
%     |            \      |
%     |              d(n) |

%==============================================================================

function D = sdiag(d)
if nargin == 0,
  help(mfilename); d = [1,2,3], D = sdiag(d), FD = full(D), return;
end;


	D = diag(sparse(d(:)));