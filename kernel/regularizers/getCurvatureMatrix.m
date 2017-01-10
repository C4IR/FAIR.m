%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function B = getCurvatureMatrix(omega,m)
%
% Builds the curvature matrix B for a domain defined by omega
% and a cell-centered grid discretization defined by m:
%
%     | \Delta          0           0 | 
% B = |      0    \Delta        0 | , for dim = 3
%     |      0       0   \Delta |
%
%
% where \Delta denotes the Laplacian operator, i.e. 
% \partial_{xx} + \partial_{yy} (2D)
% and Neumann boundary conditions are used.
%
% see also curvature.m
% =============================================================================

function B = getCurvatureMatrix(omega,m)

if nargin == 0,
  help(mfilename);
  runMinimalExample;
  return;
end;

dim = length(omega)/2;
Dxx = @(j) delta(omega,m,j);     % \partial_xx operator
switch dim,
  case 2,
    D = kron(Dxx(2),speye(m(1))) + kron(speye(m(2)),Dxx(1));
  case 3
    D = kron(Dxx(3),kron(speye(m(2)),speye(m(1)))) ...
      + kron(speye(m(3)),kron(Dxx(2),speye(m(1)))) ...
      + kron(speye(m(3)),kron(speye(m(2)),Dxx(1)));
  otherwise
    error('nyi');
end;

B = kron(speye(dim),D);

function D = delta(omega,m,j)
h = (omega(2:2:end)-omega(1:2:end))./m;
%     | -2  1       |
% D = |  1 -2  .    | / h^2
%     |     .  .  1 |
%     |        1 -2 |
D = spdiags(ones(m(j),1)*[1,-2,1],-1:1,m(j),m(j))/(h(j)*h(j));
% !! does not match descriptzion, but is OK, soso
%    D(1,1) = D(1,1)/2; D(end,end) = D(end,end)/2;
%    D([1,end]) = D([2,end-1]);

%with no BC
% D = spdiags(ones(m(j),1)*[1,-2,1],-1:1,m(j),m(j))/(h(j)*h(j));
% D([1,end],:) = 0;

%Neumann BC
D([1,end]) = -D([2,end-1]);

%------------------------------------------------------------------------------
function runMinimalExample
omega = [0,1,0,2,0,3]; m = [4,5,6];
B2 = getCurvatureMatrix(omega(1:4),m(1:2));
B3 = getCurvatureMatrix(omega,m);
figure(1); clf;
subplot(1,2,1); spy(B2); title(sprintf('%s: B2',mfilename));
subplot(1,2,2); spy(B3); title(sprintf('%s: B3',mfilename));
%==============================================================================

