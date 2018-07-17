%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% 
% Tutorial for FAIR: 3D data ordering
%
%==============================================================================

Tij = reshape(1:12,[3,2,2]), 
disp(flipdim(flipdim(permute(Tij,[3,2,1]),3),1))
disp(['T(:)''= [ 9  3 12  6  8  2 11  5  7  1 10  4]'])
%==============================================================================
