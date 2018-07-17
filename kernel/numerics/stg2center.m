%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% (c) Jan Modersitzki 2011/01/04, see FAIR.2011 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% function Y = stg2center(Y,m)
%
% transfers a staggered grid (S) to a cell-centered grid (C) or vice-versa
%
% Input:
%   Y        input points,  S or C
%   m        number of discretization points/cells
%
% Output:
%   Y        output points, C or S
%
% Example:
%  omega = [0 1 0 1]; m = [3,4];
%  yC = getCellCenteredGrid(omega,m);
%  yS = getStaggeredGrid(omega,m);
%  P = stg2center(m);
%  yC-P*yS==0, P'*yC-yS==0, yC-stg2center(yS,m)==0, stg2center(yC,m)-yS==0 
%
% see also getCellCenteredGrid, getStaggeredGrid, nodal2center, center
%==============================================================================

function P = stg2center(y,m)

if nargin == 0,
  help(mfilename);
  FAIRfigure(1); clf; 
  spy(stg2center([3,4])); 
  title(mfilename);
  return;
end;

if nargin == 1, m = y;  end;

dim = length(m);
if dim == 1,
  % note: this is nodal to cell-centered
  if nargin == 1,
    P = nodal2center(m);
  else
    P = nodal2center(y,m);
  end;
  return
end;

% ms are the sizes of the staggered components, e.g. dor dim == 3:
% 
%      | m(1)+1 m(2)   m(3)   |
% ms = | m(1)   m(2)+1 m(3)   |
%      | m(1)   m(2)   m(3)+1 |
ms  = ones(dim,1)*m + eye(dim);
ns  = prod(ms,2);
n0  = 0;

if nargin == 1,

  % ----------------------------------
  % build the matrix P and return it
  % ----------------------------------

  % for the dimension ms(i) of the i-th staggered grid, the i-th component 
  % of of m=[m(1),m(2),m(3)] has to be increased by one; the i-th staggered
  % grid contains ns(i)points

  % an averaging operator a:\R^{m+1}->\R^m is defined
  %     | 1 1     |
  % a = |   . .   |/2
  %     |     1 1 |
  a = @(j) spdiags(ones(m(j),1)*[1,1],0:1,m(j),m(j)+1)/2;

  % shortcuts for the identity matrix of size m(i) and n-by-ns(i) zeros
  E = @(i) speye(m(i),m(i));
  Z = @(i) sparse(prod(m),ns(i));

  switch dim,
    case 1,
      P = a(1);
    case 2,
      P = [kron(E(2),a(1)),  Z(2)
           Z(1),  kron(a(2),E(1))];
    case 3,
      P = [kron(E(3),kron(E(2),a(1))),  Z(2),  Z(3)
           Z(1),  kron(E(3),kron(a(2),E(1))),  Z(3)
           Z(1),  Z(2),  kron(a(3),kron(E(2),E(1))) ];
  end;
  return;
end;

% ----------------------------------
% Here starts the matrix free code
% ----------------------------------

if numel(y) == dim*prod(m),

  % cell-centered to staggered

  z = zeros(sum(ns),1);
  y = reshape(y,[],length(m));
  for i=1:length(m),
    yi = reshape(y(:,i),m);
    zi = zeros(ms(i,:));
    switch i,
      case 1, zi(1:end-1,:,:) = yi; zi(2:end,:,:) = zi(2:end,:,:) + yi;
      case 2, zi(:,1:end-1,:) = yi; zi(:,2:end,:) = zi(:,2:end,:) + yi;
      case 3, zi(:,:,1:end-1) = yi; zi(:,:,2:end) = zi(:,:,2:end) + yi;
    end;
    zi = 0.5*zi;
    z(n0+(1:ns(i))) = reshape(zi,[],1); n0 = n0 + ns(i);
  end;
else

  % staggered to cell-centered

  z = zeros(prod(m),length(m));
  for i=1:length(m),
    yi = reshape(y(n0+(1:ns(i))),ms(i,:));
    n0 = n0 + ns(i);
    switch i,
      case 1, yi = (yi(1:end-1,:,:) + yi(2:end,:,:))/2;
      case 2, yi = (yi(:,1:end-1,:) + yi(:,2:end,:))/2;
      case 3, yi = (yi(:,:,1:end-1) + yi(:,:,2:end))/2;
    end;
    z(:,i) = reshape(yi,[],1);
  end;
  z = reshape(z,[],1);
end;
P = z;
%==============================================================================
