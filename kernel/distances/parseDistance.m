%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% [m,dm] = parseDistance(type,TD,Rc,omega,m,X)
% 
% enables a generic derivative test of the distance functions which
% return a Gauss Newton type output
%==============================================================================

function [m,dm] = parseDistance(type,TD,Rc,omega,m,X)
if nargin == 0,
  m = []; dm = [];
  fprintf('void call of %s\n',mfilename)
  return;
end;
[Tc,dT] = imgModel(TD,omega,X);
[Dc,rc,dD,dr,d2psi] = distance(Tc,Rc,omega,m);
dr = dr*dT;
dD = dD*dT;

switch type,
  case 'function',      m = rc(:); dm = dr;
  case 'derivative',    m = Dc;    dm = dD;
  case 'Gauss-Newton',  m = dD';   dm = dr'*(d2psi*dr);
end;
%==============================================================================
