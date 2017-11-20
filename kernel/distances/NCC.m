%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function [Dc,rc,dD,dr,d2psi] = NCC(Tc,Rc,omega,m,varargin);
% 
% Normalized Cross Correlation based distance measure 
% using a general Gauss-Newton type framework
% computes D(Tc,Rc) =  1 - (Rc'*Tc(Y))^2/T2, Rc = Rc/norm(Rc); T2 = Tc(Y)'*Tc(Y); 
% or r = Tc, psi = 1-(r'*Rc)^2/(r'*r)
% and derivatives, dr = dT, d2psi = 2/(r'*r)
%
% Example: run('MIcc')
%   setup2DhandData;
%   xc = getCellCenteredGrid(omega,m);
%   Tc = linearInter(dataT,omega,xc);
%   Rc = linearInter(dataR,omega,xc);
%   Dc = NCC(Tc,Rc,omega,m);
%
% Input: 
%  Tc, Rc            template and reference
%  omega, m     represents domain and discretization
%  varargin        optional parameters, e.g. doDerivative=0
%
% Output:
%  Dc           NCC(Tc,Rc)
%  rc           Tc-Rc
%  dD           dpsi*dr
%  dr           dT
%  d2psi        2/(Tc'*Tc);
%
% see also distances/contents
%==============================================================================

function [Dc,rc,dD,dr,d2psi] = NCC(Tc,Rc,omega,m,varargin);

if nargin == 0,
  help(mfilename);
  setup2DhandData
  xc = getCellCenteredGrid(omega,m);
  Tc = linearInter(dataT,omega,xc);
  Rc = linearInter(dataR,omega,xc);
  D0 = feval(mfilename,Tc,Rc,omega,m);
  fprintf('%s: distance = %s\n',mfilename,num2str(D0));
  Dc = 'endOfMinimalExample'; 
  return;
end;

Dc  = []; dD = []; rc  = []; dr = []; d2psi = [];
doDerivative = (nargout > 3);

for k=1:1:length(varargin)/2, % overwrite default parameter
  eval([varargin{2*k-1},'=varargin{',int2str(2*k),'};']);
end;

%hd  = prod(omega./m); not needed here
rc = Tc; Rc = Rc/norm(Rc); T2 = Tc'*Tc;
Dc = 1-(Rc'*Tc)^2/T2;
if ~doDerivative, return; end;
dr = 1; 
dD = (-2*(Rc'*Tc)/T2*Rc'+2*((Rc'*Tc)/T2)^2*Tc')*dr;
d2psi  = 2/T2;
%==============================================================================
