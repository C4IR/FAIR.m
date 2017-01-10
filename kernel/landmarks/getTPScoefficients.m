%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function c = getTPScoefficients(LM,varargin);
%
% Computes coefficients for Thin-Plate-Spline Interpolation
% | A + theta*I   B | * | c1 |  =  | t |,
% | B'            0 |   | c2 |     | 0 |
%
% A = [rho(t_j-r_k)], B = [1,r], t=LM(:,1:dim), r = M(:,dim+1:2*dim)
%
% Input:
%   LM            landmark positions
%   varargin    optional parameter, see below
%    
% Output:
%  c            coefficients
%
% see also E5_Hands_TPS for an example
% version 2015/05/20
% =========================================================================

function c = getTPScoefficients(LM,varargin)

if nargin==0
    help(mfilename)
    runMinimalExample; 
    return;
end
theta = 0;
for k=1:2:length(varargin), % overwrites default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

% get number of landmarks, dimension, and radial basis function
n = size(LM,1);
d = size(LM,2)/2;
switch d,
  case 2,  rho   = @(r) r.^2.*log(r);
  case 3,  rho   = @(r) r;
end;
% fill 1-1 block
A = zeros(n);
for j=1:n, for k=1:j-1,
    A(j,k) = rho(norm(LM(j,d+1:end)-LM(k,d+1:end)));
end;       end;
A = A+A';
B = [ones(n,1),LM(:,d+1:end)];

% build KKT system and solve it
K = [A+theta*eye(n),B;B',zeros(d+1)];
c = K\[LM(:,1:d);zeros(d+1,d)];

%------------------------------------------------------------------------------

function runMinimalExample
setup2DhandData
c = getTPScoefficients(LM);
FAIRfigure; plot(c(:,1),'--r'); hold on; plot(c(:,2)); hold off;
title('Coefficients for Thin-Plate-Spline Interpolation')

%==============================================================================
