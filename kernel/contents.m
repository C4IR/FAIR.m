%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% This is the kernel of FAIR
%
% the kernel contains
%   contents.m	this file
%
%   data                data and examples included in the toolbox
%   distances           various distance measures (SSD etc.) and their administration
%   examples            numerous examples (organized by chapters of the FAIR book)
%   imgModel            various image models (interpolation etc.) and their administration
%   landmarks           tools for landmark registration
%   matrixfree	        tools for matrixfree operations
%   numerics            numerical tools including optimization
%   regularizers        various regularizer (elastic etc.) and their administration
%   tools               numerous tools
%   transformations     various transformation models (affine etc.) and their administration
%==============================================================================
function debit = contents
if nargout == 0, help(mfilename); return; end;

debit = {
  'contents.m'
};
%==============================================================================

