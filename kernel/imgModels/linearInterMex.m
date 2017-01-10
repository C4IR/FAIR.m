%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function [Tc,dT] = linearInterMex(T,omega,x,varargin)
%
% This wrapper runs a CPP file which codes the same scheme as linearInter
% but is hopefully faster, see linearInter for details
%==============================================================================

function [Tc,dT] = linearInterMex(T,omega,x,varargin)

Tc = mfilename('fullpath'); dT = [];
if nargin==0
    help(mfilename)
    testOneImgModel(mfilename);
    return;
elseif nargin == 1 && isempty(T)
    return;
end

% flag for computing the derivative
doDerivative = (nargout>1);
matrixFree   = 0;
boundary     = 0; % boundary condition (0: zero padding; 1: replicate)

for k=1:2:length(varargin) % overwrite default parameter
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end

dim = length(omega)/2;
m   = size(T);         if dim == 1, m = numel(T); end;
n   = length(x)/dim;

% call CPP subroutine
try
    [Tc,dT] = linearInterMexC(double(T(:)),omega,m,x(:),doDerivative,boundary==1);
catch err
    FAIRerror(err);
end

if doDerivative && not(matrixFree),
    dT = spdiags(dT,n*(0:(dim-1)),n,dim*n);
end
%==============================================================================
