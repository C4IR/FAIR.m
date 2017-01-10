%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% varargout = plotStaggeredGrid(yS,omega,m)
% 
% plots a staggered grid, returns handles to grid
% see also getStaggeredGrid, plotGrid
%==============================================================================
function varargout = plotStaggeredGrid(yS,omega,m)

if nargin == 0,
  help(mfilename)
  runMinimalExample;
  return;
end;

dim = length(omega)/2;
e  = @(i) (1:dim == i); % i-th unit vector
ns = cumsum([0;prod(ones(dim,1)*m+eye(dim),2)]);

yN = reshape(grid2grid(yS,m,'staggered','nodal'),[m+1,2]);

ph = plotGrid(yN,omega,m);
hold on;

yS11 = reshape(yS(ns(1)+1:ns(2)),m+e(1));
yS12 = (yN(:,1:end-1,2)+yN(:,2:end,2))/2;
yS21 = (yN(1:end-1,:,1)+yN(2:end,:,1))/2;
yS22 = reshape(yS(ns(2)+1:ns(3)),m+e(2));
qh = plot(yS11,yS12,'m>',yS21,yS22,'m^');


if nargout> 0,
  varargout{1} = ph;
  varargout{2} = qh;
end;
%------------------------------------------------------------------------------
function runMinimalExample

omega = [0 2 0 1];
m     = [17,13];
yS = getStaggeredGrid(omega,m);
FAIRfigure(1); clf;
plotStaggeredGrid(yS,omega,m);
%==============================================================================

