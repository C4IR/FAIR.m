%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function [y,dy,d2y] = spline1D(j,x);
%
% This function computes the values of a cubic B-spline located at j for
% the argument x as well as the first and second derivative.
%==============================================================================

function [y,dy,d2y] = spline1D(j,x)


if nargin == 0, 
    help(mfilename)
	E3_bsplines; 
	return; 
	y = 'endOfMinimalExample'; 
end;

% shift x, initialize values, derivative and second derivative
x = x-j;    y = 0*x;    dy = 0*x;   d2y = 0*x;

J = find(x>=-2 & x<-1);
y(J)   = (x(J)+2).^3;
dy(J)  = 3*(x(J)+2).^2;
d2y(J) = 6*(x(J)+2);

J = find(x>=-1 & x< 0);
y(J)   = -x(J).^3 - 2*(x(J)+1).^3+6*(x(J)+1);
dy(J)  = -3*x(J).^2 - 6*(x(J)+1).^2+6;
d2y(J) = -6*x(J) - 12*(x(J)+1);

J = find(x>= 0 & x< 1);
y(J)   = -2*(1-x(J)).^3 + x(J).^3 - 6*x(J) + 6;
dy(J)  = 6*(1-x(J)).^2 + 3*x(J).^2 - 6;
d2y(J) = -12*(1-x(J)) + 6*x(J);

J = find(x>= 1 & x< 2);
y(J)   = (2-x(J)).^3;
dy(J)  = -3*(2-x(J)).^2;
d2y(J) = 6*(2-x(J));

y = y/6; dy = dy/6; d2y = d2y/6;
%==============================================================================
