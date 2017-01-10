%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% This is a testing environment for the files in the folder
% kernel/examples
% 1.
%   Based on data/contents, a list of required files is generated and it is
%   varified, that all files are present; additional files are listed.
% 2.
%   All c-files are compiled.
% 3.
%   All files are executed.
% 4.
%   Check administrative modul distance.m
% 5.
%   test explicitly several specific distances
%==============================================================================

folder = fileparts(which(mfilename));

files = testStart(folder);
testEnd
%==============================================================================

