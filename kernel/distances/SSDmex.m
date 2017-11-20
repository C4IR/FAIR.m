%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function [Dc,rc,dD,dr,d2psi] = SSDmex(Tc,Rc,omega,m,varargin)
%
% Wrapper for C file for Sum of Squared Differences based distance measure,
% see also SSD for details 
% see also distances/contents
%==============================================================================

function [Dc,rc,dD,dr,d2psi] = SSDmex(Tc,Rc,omega,m,varargin)

if nargin == 0
    help(mfilename);
    setup2DhandData
    xc = getCellCenteredGrid(omega,m);
    Tc = linearInter(dataT,omega,xc);
    Rc = linearInter(dataR,omega,xc);
    D0 = feval(mfilename,Tc,Rc,omega,m);
    fprintf('%s: distance = %s\n',mfilename,num2str(D0))
    return
end

dD = []; dr = []; d2psi = [];
doDerivative = (nargout > 3);

for k=1:2:length(varargin)      % overwrite default parameter
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end

hd = prod((omega(2:2:end)-omega(1:2:end))./m); % voxel size for integration

if ~doDerivative
    try
        [Dc,rc] = SSDmexC(Tc,Rc,hd);
    catch err
        FAIRerror(err);
    end
else
    try
        [Dc,rc,dD,dr,d2psi] = SSDmexC(Tc,Rc,hd);
    catch err
        FAIRerror(err);
    end
end

%==============================================================================
