%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: Examples for grid generation
%
%==============================================================================

omega = [0,6,0,4,0,8]
m     = [3,2,2]
xc    = getCellCenteredGrid(omega(1:2),m(1));        xc = reshape(xc,1,[])
xc    = getCellCenteredGrid(omega(1:4),m(1:2));      xc = reshape(xc,[m(1:2),2])
xc    = getCellCenteredGrid(omega(1:6),m(1:3));      xc = reshape(xc,[m(1:3),3])
%==============================================================================
