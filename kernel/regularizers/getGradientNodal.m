%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function B = getGradientNodal(omega,m);
% 
% Builds the gradient operator B for a domain defined by omega
% and a nodal grid discretization defined by m:
%
% 
% EXAMPLE FOR 2D:
%
% o--x--o--x--o--x--o    
% |     |     |     |   
% +     +     +     +  
% |     |     |     |    
% o--x--o--x--o--x--o    
% |     |     |     |       WHERE:
% +     +     +     +       'o'  - uc       (nodal)
% |     |     |     |       'x'  - d_1 uc      (staggered-2)
% o--x--o--x--o--x--o       '+'  - d_2 uc      (staggered-1)
%
% Note: the matrices are simple, it is size that matters.  
% =============================================================================

function B = getGradientNodal(omega,m)
if nargin == 0,
  help(mfilename);
  runMinimalExample;
  return;
end;

h      =  (omega(2:2:end)-omega(1:2:end))./m;
dim    =  size(omega,2)/2;
id     =  @(i) speye(m(i)+1);   % (m(i)+1) identity matrix
                             % short difference operator
dx    =  @(i) spdiags(ones(m(i),1)*[-1,1],[0,1],m(i),m(i)+1)/h(i); 
switch dim 
    case 2
      B = [kron(id(2),dx(1));kron(dx(2),id(1))];
      z = sparse(size(B,1),size(B,2));
      B = sparse([B z; z B]);
    case 3
        B = [
        kron(id(3),kron(id(2),dx(1)))
        kron(id(3),kron(dx(2),id(1)))
        kron(dx(3),kron(id(2),id(1)))
        ];
        z = sparse(size(B,1),size(B,2));
        B = sparse([B z z; z B z; z z B]);
    otherwise
        error('Dimension must be either 2 or 3.')
end
%------------------------------------------------------------------------------
function runMinimalExample

% 2D
omega = [0,3,0,5]; m = [10,13];
B2 = getGradientNodal(omega,m);
figure(1); clf;
subplot(1,2,1)
spy(B2); 
title(sprintf('%s matrix (2D)',mfilename));


% 3D
omega = [0,1,0,2,0,3]; m = [4,5,6];
B2 = getGradientNodal(omega,m);
subplot(1,2,2);
spy(B2);   
title(sprintf('%s matrix (3D)',mfilename));
%==============================================================================
