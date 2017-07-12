%==============================================================================
% This code is part of the Matlab-based toolbox 
% LagLDDDM - A Lagrangian Gauss--Newton--Krylov Solver for Mass- and 
%                        Intensity-Preserving Diffeomorphic Image Registration
% 
% For details and license info see 
% - https://github.com/C4IR/FAIR.m/tree/master/add-ons/LagLDDMM
%
% function eigLap = eigLaplacian(omega,m)
%
% builds eigenvalues of 2D/3D Laplacian
%
% Input:
%  omega   - description of computational domain
%  m       - number of cells
%
% Output:
%  eigLap  - eigenvalues (unsorted)
%
% =========================================================================

function eigLap = eigLaplacian(omega,m)

if nargin==0
    help(mfilename);
    return;
end

dim = numel(omega)/2;


% get eigenvalues of 1D operators
eigx = eigLap1D(omega(1:2),m(1));
eigy = eigLap1D(omega(3:4),m(2));

switch dim
    case 2
        eigLap = kron(ones(m(2),1),eigx)+ kron(eigy,ones(m(1),1));
    case 3
        eigz = eigLap1D(omega(5:6),m(3));
        eigLap = kron(ones(m(3)*m(2),1),eigx) ...
                 + kron(ones(m(3),1),kron(eigy,ones(m(1),1))) ...
                 + kron(eigz,ones(prod(m(1:2)),1));
end

function eigL = eigLap1D( omega, m )
dim = numel(omega)/2;
if dim~=1, error('nyi'); end

h   = abs((omega(2)-omega(1))./m);
dx = spdiags(ones(m,1)*[-1 1],[0 1],m-1,m);
d2x = dx'*dx;

% cosine transform preparation
Cx = dct(eye(m,1));
CLx = dct(full(d2x(:,1)));

llx = CLx./Cx;
eigL = llx/h.^2;
