%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Landmark Based Thin-Plate-Spline Transformation
%
% function [Y,Z] = evalTPS(LM,c,X);
%
% evaluates the Thin-Plate-Spline parameterized by LM and c 
% at points X and the landmarks Z of the reference image.
%
% See also E5_Hands_TPS for an example
%==============================================================================

function [Y,Z] = evalTPS(LM,c,X)

if nargin==0
    help(mfilename)
    runMinimalExample; 
    return;
end

dim = size(LM,2)/2;

% prepare radial basis function
switch dim,
  case 2,  rho = @(r) r.^2.*log(r+(r==0));
  case 3,  rho = @(r) r;
end;

% L is number of landmarks, reshape X and allocate memory for 
% Y=TPS(X) and Z = TPS(r)

L = size(LM,1); 
K = dim+1:2*dim; 
Z = zeros(size(LM,1),dim);
X = reshape(X,[],dim); 
x = @(i) X(:,i); 
Y = zeros(size(X)); 

% compute the linear part
for i=1:dim,
  Y(:,i) = c(L+1,i)+X*c(L+2:end,i);
  Z(:,i) = c(L+1,i)+LM(:,K)*c(L+2:end,i);
end;

% shortcut for ||x-r_j|| in vectorized version
r = @(X,j) sqrt(sum((X-ones(size(X))*diag(LM(j,K))).^2,2));

% add non-linear part
for j=1:L,
  rj = r(X,j); Rj = r(LM(:,K),j);
  for i=1:dim,
    Y(:,i) = Y(:,i) + c(j,i)*rho(rj);   Z(:,i) = Z(:,i) + c(j,i)*rho(Rj);
  end;    
end;
Y = Y(:);Z = Z(:);

%------------------------------------------------------------------------------

function runMinimalExample
E5_2D_TPS
%==============================================================================


