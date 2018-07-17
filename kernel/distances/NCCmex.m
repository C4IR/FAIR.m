%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function [Dc,rc,dD,dr,d2psi] = NCCmex(Tc,Rc,omega,m,varargin)
%
% Wrapper for C file for Normalized Cross Correlation based distance measure,
% see also NCC for details 
% see also distances/contents
%==============================================================================

function [Dc,rc,dD,dr,d2psi] = NCCmex(Tc,Rc,omega,m,varargin)

if nargin == 0
    help(mfilename);
    setup2DhandData
    xc = getCellCenteredGrid(omega,m);
    Tc = linearInter(dataT,omega,xc);
    Rc = linearInter(dataR,omega,xc);
    D0 = feval(mfilename,Tc,Rc,omega,m);
    fprintf('%s: distance = %s\n',mfilename,num2str(D0))
    Dc = 'endOfMinimalExample';     
    return
end

dD = []; dr = []; d2psi = [];
doDerivative = (nargout > 3);

for k=1:1:length(varargin)/2  % overwrite default parameter
    eval([varargin{2*k-1},'=varargin{',int2str(2*k),'};']);
end

if ~doDerivative
    try
    [Dc,rc] = NCCmexC(Tc,Rc);
    catch err
        FAIRerror(err);
    end
else
    try
    [Dc,rc,dD,dr,d2psi] = NCCmexC(Tc,Rc);
    catch err
        FAIRerror(err);
    end       
end
%==============================================================================
