%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% function x = getStaggeredGrid(omega,m)
%
% creates a staggered grid discretization of 
%    Omega = (omega(1),omega(2)) x ... x (omega(2*dim-1,omega(2*dim)), 
%    spatial dimension dim = length(omega/2)
%
%
%          X1                     X2
% +-----+-----+-----+    +--x--+--x--+--x--+
% |     |     |     |    |     |     |     |
% x     x     x     x    |     |     |     |
% |     |     |     |    |     |     |     |
% +-----+-----+-----+    +--x--+--x--+--x--+
% |     |     |     |    |     |     |     |
% x     x     x     x    |     |     |     |
% |     |     |     |    |     |     |     |
% +-----+-----+-----+    +--x--+--x--+--x--+
%                   
% Input:
%    omega    domain specification
%   m       number of discretization points
%
% Output:
%   x       coordinates of grid points, 
%           x     = getStaggeredGrid(omega,m);
%             xS1 = reshape(xS(1:(m(1)+1)*m(2)),m(1)+1,m(2));
%             yS2 = reshape(xS((m(1)+1)*m(2)+1:end),m(1),m(2)+1);
%
% See also getCellCenteredGrid, getNodalGrid, E3_getCellCenteredGrid
%==============================================================================


function X = getStaggeredGrid(omega,m,dir)

if nargin == 0, 
  help(mfilename); 
  runMinimalExample;
  return; 
end;

% -----------------------------------------------------------------------------

X  = []; x1 = []; x2 = []; x3 = [];
h   = (omega(2:2:end)-omega(1:2:end))./m;                % voxel dimensions
% cell centers
xi = @(i) linspace(omega(2*i-1)+h(i)/2,omega(2*i)-h(i)/2,m(i))'; 
% cell boundaries (nodal)
nu = @(i) linspace(omega(2*i-1),omega(2*i),m(i)+1)'; 
switch length(omega)/2;,
  case 1, x1 = nu(1);
  case 2, 
    [x1,dum] = ndgrid(nu(1),xi(2));
    [dum,x2] = ndgrid(xi(1),nu(2));
  case 3, 
    [x1,dum,dum] = ndgrid(nu(1),xi(2),xi(3));
    [dum,x2,dum] = ndgrid(xi(1),nu(2),xi(3));
    [dum,dum,x3] = ndgrid(xi(1),xi(2),nu(3));
end;
X = [x1(:);x2(:);x3(:)];

% -----------------------------------------------------------------------------

function runMinimalExample;
omega     = [0 2 0 1]; 
m         = [6,3];
xN      = reshape(getNodalGrid(omega,m),[m+1,2]);
xS      = getStaggeredGrid(omega,m);
xS1     = reshape(xS(1:(m(1)+1)*m(2)),m(1)+1,m(2));
yS2     = xS((m(1)+1)*m(2)+1:end);
yS1     = 0.5*(xN(:,1:end-1,2)+xN(:,2:end,2));
xS2     = 0.5*(xN(1:end-1,:,1)+xN(2:end,:,1));

FAIRfigure(1,'figname',mfilename); clf;
plot(omega([1,1,2,2,1]),omega([3,4,4,3,3]),'k-','linewidth',2); hold on
plotGrid(xN,omega,m); 
plot(xS1(:),yS1(:),'r>',xS2(:),yS2(:),'r^'); axis image
title(mfilename); 

%==============================================================================
