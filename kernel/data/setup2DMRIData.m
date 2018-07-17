%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% Initializing MRI slices of a human head (T1/T2).
%==============================================================================

checkDataFile
if expfileExists, return; end;

expfile = jpgs2data('','MRIhead-T.jpg','MRIhead-R.jpg',...
  'omegaT',[0,20,0,20],'omegaR',[0,20,0,20]);
checkDataFile
%==============================================================================
