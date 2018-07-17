%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function [yc,y1,y2,y3] = grid2grid(yc,m,in,out)
%
% transfers input grid (c,s,n) to output grid (c,s,n) using inter- and extrapolation
% central operation are reduce/expand, acting on each dimension individually
%
% Input:
%   yc        input grid
%   m         number of discretization points
%   in        type of input grid
%   out       type of output grid
% Output:
%   yc        output grid
%   y1,y2,y3  components of yc
%
%==============================================================================

function [yc,y1,y2,y3] = grid2grid(yc,m,in,out)

if nargin == 0, % help and minimal example
  runMinmalExample;
  return;
end;

dim = length(m);  y{1} = []; y{3} = []; y{3} = [];
if ~exist('out','var'), out = in;  end;

switch in,
  case {'centered','cell-centered'}, n = prod(m);
    % decompose Yin
    for j=1:dim, y{j} = reshape(yc((j-1)*n+(1:n)),m); end;
    switch out,
      case {'centered','cell-centered'},
      case 'staggered', for j=1:dim, y{j} = extend(y{j},j); end;
      case 'nodal',
        for j=1:dim,
          for k=1:dim, y{j} = extend(y{j},k); end;
        end;
    end;
    
  case 'staggered',
    % decompose Yin, get dimensions right
    E = eye(dim) + m'*ones(1,dim); n = cumsum([0,prod(E)]);
    for j=1:dim, y{j} = reshape(yc(n(j)+1:n(j+1)),E(:,j)'); end;
    switch out,
      case {'centered','cell-centered'}, for j=1:dim, y{j} = reduce(y{j},j); end;
      case 'staggered',
      case 'nodal',
        for j=1:dim,
          for k=setdiff(1:dim,j), y{j} = extend(y{j},k); end;
        end;
    end;
    
  case 'nodal', n = prod(m+1);
    % decompose Yin
    for j=1:dim, y{j} = reshape(yc((j-1)*n+(1:n)),m+1); end;
    switch out,
      case {'centered','cell-centered'},
        for j=1:dim, for k=1:dim, y{j} = reduce(y{j},k); end; end;
      case 'staggered',
        for j=1:dim, for k=setdiff(1:dim,j) y{j} = reduce(y{j},k); end; end;
      case 'nodal',
    end;
    
  otherwise,
    error('can not interpret flags')
end;
y1 = y{1}; y2 = y{2}; y3 = y{3}; yc = [y1(:);y2(:);y3(:)];

%------------------------------------------------------------------------------

% the following operators act on x=y{k} in the j-th direction,
% permuting x makes j the first direction, averaging (extend creates two
% additional outside points based on linear BC) and permute back

function x = reduce(x,j)
m = size(x); J = [j,setdiff(1:length(size(x)),j)];
x = permute(x,J); x = (x(1:end-1,:,:)+x(2:end,:,:))/2; x = ipermute(x,J);

function x = extend(x,j)
m = size(x); J = [j,setdiff(1:length(size(x)),j)]; x = permute(x,J);
x = [1.5*x(1,:,:)-0.5*x(2,:,:);(x(1:end-1,:,:)+x(2:end,:,:))/2;...
  1.5*x(end,:,:)-0.5*x(end-1,:,:)];
x = ipermute(x,J);

%------------------------------------------------------------------------------

function runMinmalExample
help(mfilename);
omega = [0,6,0,4]; m = [4 3];

h   = (omega(2:2:end)-omega(1:2:end))./m;
eps = 0.2*min(h);

yC  = getCellCenteredGrid(omega,m); %yC = yC + eps*randn(size(yC));
yN  = reshape(getNodalGrid(omega,m),[],2);
yS  = getStaggeredGrid(omega,m);

yCN = grid2grid(yC,m,'centered','nodal');
yCS = grid2grid(yC,m,'centered','staggered');
ySN = grid2grid(yS,m,'staggered','nodal');
ySC = grid2grid(yS,m,'staggered','centered');
yNC = grid2grid(yN,m,'nodal','centered');
yNS = grid2grid(yN,m,'nodal','staggered');

yC  = reshape(yC, [],2);
yN  = reshape(yN, [],2);
yCN = reshape(yCN,[],2);
ySN = reshape(ySN,[],2);
ySC = reshape(ySC,[],2);
yNC = reshape(yNC,[],2);


FAIRfigure(1); clf;

subplot(3,2,1);
plotGrid(yN,omega,m,'color','g'); hold on;
plot(yC(:,1),yC(:,2),'bs');
plot(yCN(:,1),yCN(:,2),'m*');
title(sprintf('%s: example cell-centered to nodel',mfilename)); 

subplot(3,2,2);
plotGrid(yN,omega,m,'color','g'); hold on;
plot(yC(:,1),yC(:,2),'bs');
[y11,y12,y21,y22] = splittStaggeredGrid(yCS,omega,m);
plot(y11,y12,'m>',y21,y22,'m^');
title(sprintf('%s: example cell-centered to staggered',mfilename)); 

[y11,y12,y21,y22] = splittStaggeredGrid(yS,omega,m);
subplot(3,2,3);
plotGrid(yN,omega,m,'color','g'); hold on;
plot(y11,y12,'b>',y21,y22,'b^');
plot(ySC(:,1),ySC(:,2),'ms');
title(sprintf('%s: example staggered to cell-centered',mfilename)); 

subplot(3,2,4);
plotGrid(yN,omega,m,'color','g'); hold on;
plot(y11,y12,'b>',y21,y22,'b^');
plot(ySN(:,1),ySN(:,2),'m*');
title(sprintf('%s: example staggered to nodal',mfilename)); 

subplot(3,2,5);
plotGrid(yN,omega,m,'color','g'); hold on;
plot(yN(:,1),yN(:,2),'b*');
plot(yNC(:,1),yNC(:,2),'ms');
title(sprintf('%s: example nodal to cell-centered',mfilename)); 

subplot(3,2,6);
plotGrid(yN,omega,m,'color','g'); hold on;
plot(yN(:,1),yN(:,2),'b*');
[y11,y12,y21,y22] = splittStaggeredGrid(yNS,omega,m);
plot(y11,y12,'m>',y21,y22,'m^');

title(sprintf('%s: example nodal to staggered',mfilename)); 

%------------------------------------------------------------------------------

function [y11,y12,y21,y22] = splittStaggeredGrid(yS,omega,m);
dim = length(omega)/2;
e  = @(i) (1:dim == i); % i-th unit vector
ns = cumsum([0;prod(ones(dim,1)*m+eye(dim),2)]);
yN = reshape(grid2grid(yS,m,'staggered','nodal'),[m+1,2]);
y11 = reshape(yS(ns(1)+1:ns(2)),m+e(1));
y12 = (yN(:,1:end-1,2)+yN(:,2:end,2))/2;
y21 = (yN(1:end-1,:,1)+yN(2:end,:,1))/2;
y22 = reshape(yS(ns(2)+1:ns(3)),m+e(2));

%==============================================================================
