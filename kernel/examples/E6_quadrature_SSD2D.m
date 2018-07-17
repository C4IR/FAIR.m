%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR:  SSD versus discretization width (midpoint quadrature rule)
%
% - load data (setup2DHNSPData)
% - interpolate on various grids (splineInter)
% - compute SSD and plot it
%==============================================================================

clear; help(mfilename);

setup2DHNSPData; clf; h = []; Q = [];
imgModel('reset','imgModel','splineInter');
[T,R] = imgModel('coefficients',dataT,dataR,omega,'out',0);
for j=1:10,
  m    = 2^j*[1,1]; 
  h(j) = prod((omega(2:2:end)-omega(1:2:end))./m); 
  xc   = getCellCenteredGrid(omega,m); 
  res  = imgModel(T,omega,xc) - imgModel(R,omega,xc);
  psi  = 0.5*h(j)*res'*res;
  Q(j) = psi;
end;
figure(1); clf; p1=semilogx(h/h(1),Q+eps,'kx',h/h(1),Q,'k-');
%==============================================================================
