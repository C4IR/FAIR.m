%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% function X = getNodalGrid(omega,m)
% (c) Jan Modersitzki 2009/04/08, see FAIR.2 and FAIRcopyright.m.
%
%
% creates a nodal discretization of 
%    Omega = (omega(1),omega(2)) x ... x (omega(2*dim-1,omega(2*dim)), 
%    spatial dimension dim = length(omega/2)
%
%
%          X1                     X2
% x-----x-----x-----x    x-----x-----x-----x
% |     |     |     |    |     |     |     |
% |     |     |     |    |     |     |     |
% |     |     |     |    |     |     |     |
% x-----x-----x-----x    x-----x-----x-----x
% |     |     |     |    |     |     |     |
% |     |     |     |    |     |     |     |
% |     |     |     |    |     |     |     |
% x-----x-----x-----x    x-----x-----x-----x
%
% Input:
%    omega    domain specification
%   m       number of discretization points
%
% Output:
%   x       coordinates of grid points, 
%           x is d*prod(m+1)-by-1
%                              | x^1_1, ... x^dim_1 |
%                              |        ...         |
%            reshape(x,dim) is | x^1_j, ... x^dim_j |
%                              |        ...         |
%                              | x^1_n, ... x^dim_n |
%             with x_j = (x^1_1,...,x^dim_j) is the j-th grid point, 
%            j = 1,...,n = prod(m+1), dim = length(omega)/2
%
% See also getCellCenteredGrid, getStaggeredGrid, E3_getCellCenteredGrid
%==============================================================================

function X = getNodalGrid(omega,m)

if nargin == 0, % help and minimal example
  help(mfilename); 
  runMinimalExample;
  return; 
end;

% -----------------------------------------------------------------------------

X  = []; x1 = []; x2 = []; x3 = [];
h   = (omega(2:2:end)-omega(1:2:end))./m; % voxel size for integration
% cell boundaries (nodal)
nu = @(i) linspace(omega(2*i-1),omega(2*i),m(i)+1)';   
switch length(omega)/2,
  case 1, x1 = nu(1);
  case 2, [x1,x2] = ndgrid(nu(1),nu(2));
  case 3, [x1,x2,x3] = ndgrid(nu(1),nu(2),nu(3));
end;
X = [x1(:);x2(:);x3(:)];

% -----------------------------------------------------------------------------

function runMinimalExample;
  
omega     = [0 2 0 1]; 
m         = [6,3];
xN         = getNodalGrid(omega,m);
xP         = reshape(xN,[],2)';

FAIRfigure(1,'figname',mfilename); clf;
plot(omega([1,1,2,2,1]),omega([3,4,4,3,3]),'k-','linewidth',2); hold on
plotGrid(xN,omega,m); 
plot(xP(1,:),xP(2,:),'rs'); axis image
title(mfilename); 

%==============================================================================
