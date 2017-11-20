%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
% ##2
%==============================================================================
%
% function [Tc,dT] = linearInterSmoothMex(T,omega,x,varargin);
%
% This wrapper runs a CPP file which codes the same scheme as linearInterSmooth 
% but is hopefully faster, see linearInterSmooth for details
%==============================================================================

function [Tc,dT] = linearInterSmoothMex(T,omega,x,varargin)

Tc = mfilename('fullpath'); dT = [];
if nargin == 0
    testOneImgModel(mfilename);
    return
elseif nargin == 1 && isempty(T)
    return
end

eta = 0.1;
% flag for computing the derivative
doDerivative = (nargout>1);
matrixFree   = 0;

for k=1:2:length(varargin) % overwrite default parameter
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end

% get data size m, dimension d, and number n of interpolation points
dim = length(omega)/2;
m   = size(T);         if dim == 1, m = numel(T); end;
n   = length(x)/dim;

% call CPP subroutine
try
    [Tc,dT] = linearInterSmoothMexC(double(T(:)),omega,m,x(:),eta,doDerivative);
catch err
    FAIRerror(err);
end
if doDerivative && not(matrixFree),
    dT = spdiags(dT,n*(0:(dim-1)),n,dim*n); 
end
%==============================================================================
