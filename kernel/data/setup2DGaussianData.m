%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% Artificial data motivated by mass-preserving imaging
%==============================================================================

checkDataFile
if expfileExists, return; end;

viewPara = {'viewImage','viewImage2D','colormap',gray(256)};
imgPara  = {'imgModel','splineInter','regularizer','moments','theta',1e-2};
traPara  = {'trafo','affine2D'};
disPara  = {'distance','SSD'};
regPara  = {'regularizer', 'mbHyperElastic','alpha',1,...
  'alphaLength',1000,'alphaArea',0,'alphaVolume',100};

expfile = jpgs2data('','Gauss-T.jpg','Gauss-R.jpg', ...
  'omega',[0,1,0,1],'m',[128,128],...
  'viewPara',viewPara,'imgPara',imgPara);

save(expfile,'-append','viewPara','imgPara','traPara','disPara','regPara');
checkDataFile
%==============================================================================
