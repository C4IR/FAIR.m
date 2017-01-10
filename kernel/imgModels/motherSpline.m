%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function b = motherSpline(j,xi)
%
% shortcut for the mother spline function on the four non-trivial 
% intervals, j \in{1,2,3,4}
% and its derivatives, j\in{5,6,7,8}, see FAIR, p.27 Eq. (3.7).
%==============================================================================

function b = motherSpline(j,xi)

if nargin == 0,
  help(mfilename);
  runMinimalExample;
  return;
end;

switch j,
  case 1,  b = (2+xi).*(2+xi).*(2+xi);  % MATLAB! b = (2+xi).^3;
  case 2,  b = -(3*xi+6).*xi.^2+4;      % -xi.^3 - 2*(xi+1).^3 +6*(xi+1);
  case 3,  b =  (3*xi-6).*xi.^2+4;      % xi.^3 + 2*(xi-1).^3 -6*(xi-1);
  case 4,  b = (2-xi).*(2-xi).*(2-xi);  % MATLAB! b = (2-xi).^3;
  case 5,  b =  3*(2+xi).^2;
  case 6,  b = -(9*xi+12).*xi;          % -3*xi.^2 - 6*(xi+1).^2+6;
  case 7,  b =  (9*xi-12).*xi;          %  3*xi.^2 + 6*(xi-1).^2-6;
  case 8,  b = -3*(2-xi).^2;
end;

%------------------------------------------------------------------------------

function runMinimalExample

fprintf('%s: minimal example\n',mfilename)
m = 101;  
x = linspace(0,1,m);
z = zeros(m,4);  
y = zeros(m,4,2);  
for j=1:4,
  z(:,j)   = x+j-3;
  y(:,j,1) = motherSpline(j,  z(:,j));
  y(:,j,2) = motherSpline(j+4,z(:,j));
end;

FAIRfigure(1); clf; 
plot(z(:),reshape(y(:,:,1),1,[]),'b-',z(:),reshape(y(:,:,2),1,[]),'r-'); 
title(mfilename);
legend('motherspline','derivative','location','southwest')
title(mfilename);
%==============================================================================


