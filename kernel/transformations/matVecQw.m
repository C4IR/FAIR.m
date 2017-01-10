%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function Qw = matVecQw(w,Q,flag)
%
% computes y = Q*w or Q'*w
% where Q =kron(Q3,kron(Q2,Q1)) is efficiently stored as Q= {Q1,Q2,Q3}
% see also splineTransformation2Dsparse 
%==============================================================================

function Qw = matVecQw(w,Q,flag)

if nargin == 0,
  help(mfilename);
  return;
end;

% default is Q*w
if ~exist('flag','var'), flag = 'Qw';  end;

if strcmp(flag,'QTw'),
  % build transposes, note: small matrices
  for i=1:length(Q),    Q{i}=Q{i}';  end;
  Qw = matVecQw(w,Q);
  return;
end;

% reconstruct sizes m and p from Q and w, Q{j} is m{j}-by-p{j}

dim = length(Q);
for j=1:length(Q),
  p(j) = size(Q{j},2);
  m(j) = size(Q{j},1);
end;

switch dim,
  case 1,
    Qw = Q{1}*w;
  case 2,
    w  = reshape(w,[p,dim]);
    Qw = zeros(m(1),m(2),dim);
    for i =1:dim,
      Qw(:,:,i) = Q{1}*w(:,:,i)*Q{2}'; 
    end;
  case 3,
    w  = reshape(w,[p,dim]);
    Qw = zeros(m(1),m(2),m(3),dim);
    for i =1:dim,
        Qw(:,:,:,i) = tensorProd(Q,w(:,:,:,i));
    end;    
end;

Qw = reshape(Qw,[],1);
%------------------------------------------------------------------------------

function Qw = tensorProd(Q,w);
for i=1:length(Q),
  p(i,:) = size(Q{i});
end;
Qw = zeros(p(:,1)');

for i1=1:size(Qw,1),
  for i2=1:size(Qw,2),
    for i3=1:size(Qw,3),
      for j1=1:size(w,1),
        for j2=1:size(w,2),
          for j3=1:size(w,3),
            Qw(i1,i2,i3) = Q{1}(i1,j1)*Q{2}(i2,j2)*Q{3}(i3,j3)*w(j1,j2,j3);
          end;
        end;
      end;
    end;
  end;
end;
%==============================================================================

