%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% PET-CT data of a human thorax. Original data from
%     @article{ShekharEtAl2005,
%      author = {Raj Shekhar and Vivek Walimbe and Shanker Raja and Vladimir Zagrodsky
%                and Mangesh Kanvinde and Guiyun Wu and Bohdan Bybel},
%      title = {Automated 3-Dimensional Elastic Registration of Whole-Body {PET}
%               and {CT} from Separate or Combined Scanners},
%      journal = {J. of Nuclear Medicine},
%      volume = {46},
%      number = {9},
%      year = {2005},
%      pages = {1488--1496},
%     }
%==============================================================================

checkDataFile
if expfileExists, return; end;

expfile = jpgs2data('','PET-CT-PET.jpg','PET-CT-CT.jpg', ...
  'omega',[0,50,0,50],'m',[128,128]);
checkDataFile
%==============================================================================

