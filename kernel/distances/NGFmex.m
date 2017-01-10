%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function [Dc,rc,dD,dr,d2psi] = NGFmex(Tc,Rc,omega,m,varargin)
%
% Wrapper for C file for Normalized Gradient Field based distance measure,
% see also NGFdot for details 
% see also distances/contents
%==============================================================================

function [Dc,rc,dD,drc,d2psi] = NGFmex(Tc,Rc,omega,m,varargin)

if nargin == 0
   help(mfilename);
   setup2DhandData
   xc = getCellCenteredGrid(omega,m);
   Tc = linearInter(dataT,omega,xc);
   Rc = linearInter(dataR,omega,xc);
   D0 = feval(mfilename,Tc,Rc,omega,m);
   fprintf('%s: distance = %s\n',mfilename,num2str(D0));
   Dc = 'endOfMinimalExample';
   return;
end

dD = []; drc = []; d2psi = [];
edge         = 100;                       % the edge parameter
doDerivative = (nargout > 2);             % compute derivatives only on demand

for k=1:2:length(varargin)                % overwrites default parameter
   eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end

o = prod(omega(2:2:end)-omega(1:2:end));
h = (omega(2:2:end)-omega(1:2:end))./m;   % voxel size for integration

if ~doDerivative
    try
        [Dc,rc] = NGFdotMexC(Tc,Rc,o,m,h,edge);
    catch err
        FAIRerror(err);
    end
else   
   [Dc,rc,dD,drc,d2psi] = NGFdotMexC(Tc,Rc,o,m,h,edge);
end
%==============================================================================

