%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% converts a cell into a struct and vice versa
%==============================================================================

function out = FAIRcell2struct(in)

if nargin == 0,
  return;
end;

if isstruct(in),
  f = fieldnames(in);
  v = struct2cell(in);
  out = reshape({f{:};v{:}},1,[]);
elseif iscell(in)
  out = cell2struct(in(2:2:end),in(1:2:end),2);
else
  fprintf(2,'can not deal input type\n');
  whos
  error(1)
end;
%==============================================================================

  