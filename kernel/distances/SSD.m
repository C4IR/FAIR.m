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
% For use of weights see line 38. 
%
% computes D(Tc,Rc) = hd*psi(r(Tc)), r = Tc-Rc, psi = 0.5*r'*W*r 
% and derivatives, dr = dT, d2psi= hd*W, hd = prod(omega./m). 
% W either scalar or mastrix, e.g., W=speye(prod(m)); details see below. 
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
%  varargin   optional parameters, e.g. doDerivative=0, weights=1.0
%
% Output:
%  Dc         SSD(Tc,Rc)
%  rc         Tc-Rc
%  dD         dpsi*dr
%  dr         dT
%  d2psi      d2psi = prod(omega./m)*diag(weights)
% 
% WEIGHTING: W is either scalar or matrix. For an example using weighted SSD 
%            in multi-level registration see E9_Hands_MLIR_wSSD_mfElas.m
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

Dc  = []; dD = []; rc  = []; dr = []; 
doDerivative = (nargout > 2);
weights = 1.0;
for k=1:2:length(varargin), % overwrite default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

hd = prod((omega(2:2:end)-omega(1:2:end))./m); % voxel size for integration

d2psi = hd*weights;

rc = Tc - Rc;                         % the residual
Dc = 0.5 * (rc'*d2psi*rc);            % the SSD

if ~doDerivative, return; end;

dr = 1;                         % or speye(length(rc),length(rc)); 
dD = rc'*d2psi*dr; 

%------------------------------------------------------------------------------
function runMinimalExample

setup2DhandData;
xc = getCellCenteredGrid(omega,m);
Tc = linearInter(dataT,omega,xc);
Rc = linearInter(dataR,omega,xc);
D0 = feval(mfilename,Tc,Rc,omega,m,'doDerivative',1);
fprintf('%s: distance = %s\n',mfilename,num2str(D0))

%==============================================================================
