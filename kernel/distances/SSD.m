%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function [Dc,rc,dD,dr,d2psi] = SSD(Tc,Rc,omega,m,varargin);
%
% Sum of Squared Differences based distance measure 
% for usage in a general Gauss-Newton type framework
%
% computes D(Tc,Rc) = hd*psi(r(Tc)), r = Tc-Rc, psi = 0.5*r'*r 
% and derivatives, dr = dT, d2psi= hd*I, hd = prod(omega./m)
%
%   setup2DhandData;
%   xc = getCellCenteredGrid(omega,m);
%   Tc = linearInter(dataT,omega,xc);
%   Rc = linearInter(dataR,omega,xc);
%   Dc = SSD(Tc,Rc,omega,m);
%
% Input: 
%  Tc, Rc     template and reference
%  omega, m   represents domain and discretization
%  varargin   optional parameters, e.g. doDerivative=0
%
% Output:
%  Dc         SSD(Tc,Rc)
%  rc         Tc-Rc
%  dD         dpsi*dr
%  dr         dT
%  d2psi      hd = prod(omega./m)
% 
% see also distances/contents
%==============================================================================

function [Dc,rc,dD,dr,d2psi] = SSD(Tc,Rc,omega,m,varargin);

if nargin==0
    help(mfilename)
    runMinimalExample; 
     Dc = 'endOfMinimalExample';     
    return;
end

Dc  = []; dD = []; rc  = []; dr = []; d2psi = [];
doDerivative = (nargout > 2);

for k=1:2:length(varargin), % overwrite default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

hd = prod((omega(2:2:end)-omega(1:2:end))./m); % voxel size for integration
rc = Tc - Rc;                         % the residual
Dc = 0.5*hd * (rc'*rc);            % the SSD

if ~doDerivative, return; end;

dr = 1;                         % or speye(length(rc),length(rc)); 
dD = hd * rc'*dr; 
d2psi = hd;

%------------------------------------------------------------------------------
function runMinimalExample

setup2DhandData;
xc = getCellCenteredGrid(omega,m);
Tc = linearInter(dataT,omega,xc);
Rc = linearInter(dataR,omega,xc);
D0 = feval(mfilename,Tc,Rc,omega,m,'doDerivative',1);
fprintf('%s: distance = %s\n',mfilename,num2str(D0))

%==============================================================================
