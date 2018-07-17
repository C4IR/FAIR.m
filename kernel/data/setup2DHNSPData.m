%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% Images from a histologial serial sectioning,
% data courtesy Oliver Schmitt, Institute of Anatomy, University of Rostock, Germany
%==============================================================================

checkDataFile
if expfileExists, return; end;
viewPara = {'viewImage','viewImage2D','colormap','gray(256)'};
expfile = jpgs2data('','HNSP-T.jpg','HNSP-R.jpg','omega',[0,2,0,1],'viewPara',viewPara,'m',[512,256]);
checkDataFile
%==============================================================================
