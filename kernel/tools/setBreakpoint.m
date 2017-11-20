%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% adds a list of breakpoints to an mfile
% function setBreakpoint(file,list)
%==============================================================================

function setBreakpoint(file,list)
if nargin == 0,
  help(mfilename)
  return;
end;

FAIRtestStatus = getappdata(0,'FAIRtestStatus');

if isfield(FAIRtestStatus,'FAIRrun') ...
    & strcmp(FAIRtestStatus.('FAIRrun'),'on')
  fprintf('FAIRrun == ''on'' => no breakpoints set in %s\n',mfilename);
  return;
end;


for j=1:length(list),
  eval(sprintf('dbstop in %s at %d\n',file,list{j}));
end;

%==============================================================================