%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% function x = getCellCenteredGrid(omega,m)
%
% creates a cell-center discretization of 
%    Omega = (omega(1),omega(2)) x ... x (omega(2*dim-1,omega(2*dim)), 
%    spatial dimension d = length(omega/2)
%
%          x1                     x2
% +-----+-----+-----+    +-----+-----+-----+
% |     |     |     |    |     |     |     |
% |  x  |  x  |  x  |    |  x  |  x  |  x  |
% |     |     |     |    |     |     |     |
% +-----+-----+-----+    +-----+-----+-----+
% |     |     |     |    |     |     |     |
% |  x  |  x  |  x  |    |  x  |  x  |  x  |
% |     |     |     |    |     |     |     |
% +-----+-----+-----+    +-----+-----+-----+
%                   
% Input:
%    omega    domain specification
%   m       number of discretization points
%
% Output:
%   x       coordinates of grid points, 
%           x is dim*prod(m)-by-1
%                              | x^1_1, ... x^dim_1 |
%                              |        ...         |
%            reshape(x,dim) is | x^1_j, ... x^dim_j |
%                              |        ...         |
%                              | x^1_n, ... x^dim_n |
%             with x_j = (x^1_1,...,x^dim_j) is the j-th grid point, 
%            j = 1,...,n = prod(m), dim = length(omega)/2;
% See also getStaggeredGrid, getNodalGrid, E3_getCellCenteredGrid
%==============================================================================

function x = getCellCenteredGrid(omega,m)

if nargin == 0,     % help and minimal example
  help(mfilename); 
  runMinimalExample;
  return; 
end;

% -----------------------------------------------------------------------------

x  = []; x1 = []; x2 = []; x3 = []; x4=[];
h  = (omega(2:2:end)-omega(1:2:end))./m;                             % voxel dimensions
xi = @(i) linspace(omega(2*i-1)+h(i)/2,omega(2*i)-h(i)/2,m(i))';     % cell centers

switch length(omega)/2;,
  case 1, x1 = xi(1);
  case 2, [x1,x2] = ndgrid(xi(1),xi(2));
  case 3, [x1,x2,x3] = ndgrid(xi(1),xi(2),xi(3));
  case 4, [x1,x2,x3,x4] = ndgrid(xi(1),xi(2),xi(3),xi(4));
end;
x = [reshape(x1,[],1);reshape(x2,[],1);reshape(x3,[],1); reshape(x4,[],1)];

% -----------------------------------------------------------------------------

function runMinimalExample;
omega     = [0 2 0 1]; 
m         = [6,3];
xN         = getNodalGrid(omega,m);
xC         = reshape(getCellCenteredGrid(omega,m),[],2)';

FAIRfigure(1,'figname',mfilename); clf;
plot(omega([1,1,2,2,1]),omega([3,4,4,3,3]),'k-','linewidth',2); hold on
plotGrid(xN,omega,m); 
plot(xC(1,:),xC(2,:),'rx'); axis image
title(mfilename); 

%==============================================================================
