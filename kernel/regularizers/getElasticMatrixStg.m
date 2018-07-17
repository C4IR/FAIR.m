%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function B = getElasticMatrixStg(omega,m,mu,lambda)
%
% Builds the elasticity matrix B for a domain defined by omega
% and a staggered grid discretization defined by m:
%
%     | a\nabla       0        0 | 
% B = |       0 a\nabla        0 |
%     |       0       0  a\nabla |
%     |           b \div         |,
%
% where \nabla and \div are weighted by a=\sqrt(mu) and b=\sqrt(\lambda+\mu),
% respectively; defaults: \mu = 1, \lambda = 0.
% Note: the matrices are simple, it is size that matters.
%
% =============================================================================

function B = getElasticMatrixStg(omega,m,mu,lambda)
if nargin == 0,
  help(mfilename);
  runMinimalExample;
  B = 'endOfMinimalExample';
  return;
end;

                                        % set defautls
if ~exist('mu','var'),     mu     = [];  end;
if ~exist('lambda','var'), lambda = [];  end;
if isempty(mu),     
  warning('mu has not yet been set! set mu=1');  
  mu = 1;  
end;
if isempty(lambda), 
  warning('lambda has not yet been set! lambda=0');  
  lambda = 0;  
end;

dim = length(omega)/2; 
h   = (omega(2:2:end)-omega(1:2:end))./m; % voxel size for integration
a   = sqrt(mu); 
b   = sqrt(mu+lambda);

% for the dimension ms(i) of the i-th staggered grid, the i-th component 
% of of m=[m(1),m(2),m(3)] has to be increased by one; the i-th staggered 
% grid contains ns(i) points

ms  = @(i)   m + ((1:dim)==i); 
ns  = @(i)   prod(ms(i));
% Example: m=(2,3), ms(1) = (3,3); ms(2)=(2,4), ns(1) = 9, ns(2) = 8.

% the dimension of \partial^{h,i}_k is p(i,k)-by-ns(i),
% where q(i)=prod(ms(i)) is the length of the i-th staggered grid 
% and p(i,k) is the product of a modified m: 
%   m(i) has to be increase by one and
%   m(k) has to be decrase  by one
p = @(i,k) prod(m + ((1:dim)==i) - ((1:dim==k)));

% Example:
% p(1,1) = 2*3         = 6,   p(1,2) = (2+1)*(3-1) = 6, ns(1) = 9
% p(2,1) = (2-1)*(3+1) = 4,   p(2,2) = 2*3         = 6, ns(2) = 8

% E is a q-by-q identity matrix, where q=m(k)+1 for i==k and q=m(k) (i~=k)
E = @(i,k) speye(m(k)+(i==k));

% Example: 
% E(1,1) = speye(m(1)+1) = speye(3)
% E(1,2) = speye(m(2))   = speye(3)
% E(2,1) = speye(m(1))   = speye(2)
% E(2,2) = speye(m(2)+1) = speye(4)

% D is the 1D derivative matrix of size p-by-(p+1)
%        |-1 1     |      
% D = |  .  .   |/h(k), where p=m(k) (i==k) and p=m(k)-1 (i~=k)
%        |     -1 1|      

D = @(i,k) spdiags(ones(m(k)-(i~=k),1)*[-1,1],...
  [0,1],m(k)-(i~=k),m(k)+(i==k))/h(k);

% Example: 
% D(1,1) is m(1)-by-(m(1)+1) = 2-by-3
% D(1,2) is (m(2)-1)-by-m(2) = 2-by-3
% D(2,1) is (m(1)-1)-by-m(1) = 1-by-2
% D(2,2) is m(2)-by-(m(2)+1) = 3-by-4


switch dim
  case 2
    % build the 2D elasticity operator
    %
    %      | a\nabla 0      |   |a D11  0     |
    %  B = |                | = |a D12  0     |
    %      | 0      a\nabla |   |0      a D21 |
    %      |                |   |0      a D22 |
    %      | b Div1 b Div2  |   |b D11  b D22 |

    D11 = kron(E(1,2),D(1,1));  D12 = kron(D(1,2),E(1,1));
    D21 = kron(E(2,2),D(2,1));  D22 = kron(D(2,2),E(2,1));
    
    % Example:
    % D11 = 3-by-3 \otimes 2-by-3 = 6-by-9
    % D12 = 2-by-3 \otimes 3-by-3 = 6-by-9
    % D21 = 4-by-4 \otimes 1-by-2 = 4-by-8
    % D22 = 3-by-4 \otimes 2-by-2 = 6-by-8
    
    B = [  a*D11,sparse(p(1,1),ns(2));
           a*D12,sparse(p(1,2),ns(2));
           sparse(p(2,1),ns(1)),a*D21;
           sparse(p(2,2),ns(1)),a*D22;
           b*D11,b*D22];
    % Example: B = 28-by-17
  case 3
    % build the 3D elasticity operator
    %
    %      |a D11  0      0    |
    %  B = |a D12  0      0    |
    %      |a D13  0      0    |
    %      |0      a D21  0    |
    %      |0      a D22  0    |
    %      |0      a D23  0    |
    %      |0      0      a D31|
    %      |0      0      a D32|
    %      |0      0      a D33|
    %      |b D11  b D22  b D33|

    % Example: m=[2,3,4] ns = (36,32,30)
    D11 = kron(kron(E(1,3),E(1,2)),D(1,1));
    D12 = kron(kron(E(1,3),D(1,2)),E(1,1));
    D13 = kron(kron(D(1,3),E(1,2)),E(1,1));
    D21 = kron(kron(E(2,3),E(2,2)),D(2,1));
    D22 = kron(kron(E(2,3),D(2,2)),E(2,1));
    D23 = kron(kron(D(2,3),E(2,2)),E(2,1));
    D31 = kron(kron(E(3,3),E(3,2)),D(3,1));
    D32 = kron(kron(E(3,3),D(3,2)),E(3,1));
    D33 = kron(kron(D(3,3),E(3,2)),E(3,1));   

    % Example:
    % D11 = kron(kron(4-by-4,3-by-3),2-by-3) = 24-by-36
    % D12 = kron(kron(4-by-4,2-by-3),3-by-3) = 24-by-36
    % D13 = kron(kron(3-by-4,3-by-3),3-by-3) = 27-by-36
    % D21 = kron(kron(4-by-4,4-by-4),1-by-2) = 16-by-32
    % D22 = kron(kron(4-by-4,3-by-4),2-by-2) = 24-by-32
    % D23 = kron(kron(3-by-4,4-by-4),2-by-2) = 24-by-32
    % D31 = kron(kron(5-by-5,3-by-3),1-by-2) = 15-by-30
    % D32 = kron(kron(5-by-5,2-by-3),2-by-2) = 20-by-30
    % D32 = kron(kron(4-by-5,3-by-3),2-by-2) = 24-by-30
        
    B = [
      a * D11, sparse(p(1,1),ns(2)), sparse(p(1,1),ns(3)); ...
      a * D12, sparse(p(1,2),ns(2)), sparse(p(1,2),ns(3)); ...
      a * D13, sparse(p(1,3),ns(2)), sparse(p(1,3),ns(3)); ...
      sparse(p(2,1),ns(1)), a * D21, sparse(p(2,1),ns(3)); ...
      sparse(p(2,2),ns(1)), a * D22, sparse(p(2,2),ns(3)); ...
      sparse(p(2,3),ns(1)), a * D23, sparse(p(2,3),ns(3)); ...
      sparse(p(3,1),ns(1)), sparse(p(3,1),ns(2)), a * D31; ...
      sparse(p(3,2),ns(1)), sparse(p(3,2),ns(2)), a * D32; ...
      sparse(p(3,3),ns(1)), sparse(p(3,3),ns(2)), a * D33; ...
      b * D11,              b * D22,              b * D33      
      ];
    % Example:
    % B = [75-by-36    0      0
    %         0     64-by-32  0
    %         0        0     59-by-30
    %      24-by-36 24-by-32 24-by-30 ] 
    %   = 228-by-98
    
  otherwise,  error('nyi');
end;
%------------------------------------------------------------------------------
function runMinimalExample

omega = [0,1,0,2,0,3]; m = [4,5,6];
B2 = getElasticMatrixStg(omega(1:4),m(1:2),1,0);
B3 = getElasticMatrixStg(omega,m,1,0);
figure(1); clf;
subplot(1,2,1); spy(B2); title(sprintf('%s: B2',mfilename));
subplot(1,2,2); spy(B3); title(sprintf('%s: B3',mfilename));
%==============================================================================



