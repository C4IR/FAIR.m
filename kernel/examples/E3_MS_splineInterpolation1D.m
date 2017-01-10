%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: multi-scale spline model for 1D data
%
%==============================================================================

function varargout = E3_MS_splineInterpolation1D(regularizer)
if nargin == 0, regularizer = 'I';  end;                  % default value
dataT  = [0,2,2,2,1]'; m = length(dataT); omega = [0,m];  % initialize data
xc = linspace(omega(1)-1,omega(2)+1,101);                 % fine discretization 
Tc = @(T) splineInter(T,omega,xc);                        % spline interpolant
B  = spdiags(ones(m,1)*[1,4,1],[-1:1],m,m);               % spline basis
D  = spdiags(ones(m,1)*[-1,1],[0,1],m-1,m);               % derivative operation
M  = toeplitz([96,-54,0,6,zeros(1,m-4)]);                 % second derivative 
switch regularizer,                                       % initialize regularization   
  case 'I', W  = speye(m,m);                              % identity,   Tikhonov 
  case 'D', W  = D'*D;                                    % derivative, Tikhonov-Phillips
  case 'M', W  = M;                                       % bending operator
end;
c  = @(theta) (B'*B+theta*W)\(B'*dataT);  % coefficients as function in  theta
FAIRfigure(5); clf;
ph(1:2) = plot(getCellCenteredGrid(omega,m),dataT,'.k',xc,Tc(B\dataT),'k-');  hold on;
ph(3:5) = plot(xc,Tc(c(1)),'k-.',xc,Tc(c(19)),'k--',xc,Tc(c(100)),'k-');  hold off;
if nargout ~= 0; varargout = {ph}; end;
%==============================================================================
