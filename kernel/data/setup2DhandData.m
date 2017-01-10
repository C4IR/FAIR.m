%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% Initializes xray data of human hands and adds manually chosen landmarks.
%
% Data originates from
%  @article{Amit1994,
%    author = {Yali Amit},
%     title = {A nonlinear variational problem for image matching},
%      year = {1994},
%   journal = {SIAM J. Sci. Comput.},
%    volume = {15},
%    number = {1},
%     pages = {207--224},
%  }
%==============================================================================

checkDataFile;        if expfileExists, return; end;

% for this data landmarks have been manually identified by JM
LM    = [
  5.5841   17.2664    2.6807   12.7797
  10.7243   21.6121    7.2028   19.6795
  13.2477   21.6121   10.1865   20.8916
  15.2570   19.2290   12.5175   20.0991
  15.8645   15.1636   14.3357   16.7424
  5.3972    8.1075    7.9953    6.3462
  7.5000    5.9579   11.8648    5.6469
  ];

% set special view options
viewPara = {'viewImage','viewImage2D','colormap','bone(256)'};
expfile = jpgs2data('','hands-T.jpg','hands-R.jpg',...
  'omega',[0,20,0,25],'m',[128,128],'viewPara',viewPara,'LM',LM);
checkDataFile
%==============================================================================

