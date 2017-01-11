%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% function [Dc,rho,dD,drho,d2psi] = MImex(Tc,Rc,omega,m,varargin);
%
% Mutual Information based distance measure  (using a C implemetation for the joint
% density, see MIspline for pure MATLAB)  using a general Gauss-Newton type framework.
%
% The MI distance measure is implemented using a Parzen-window estimator
% for the joint density rho. Given rho and drho, MI is computed as
% Dc    = psi(rho) = rho' * log(rho + tol) + ...
% dpsi  = log(rho + tol) + rho./(rho + tol) + ...
% d2psi = (rho + 2*tol)./(rho + tol)^2 + ...
% where the "..." indicates two more summands of similar type
%
% Example: run('MIcc')
%   setup2DHandData;
%   xc = getCellCenteredGrid(omega,m);
%   Tc = linearInter(dataT,omega,xc);
%   Rc = linearInter(dataR,omega,xc);
%   Dc = MIcc(Tc,Rc,omega,m);
%
% Input:
%  T, R            template and reference
%  omega, m     represents domain and discretization
%  varargin        optional parameters, e.g. doDerivative=0
%
% Output:
%  Dc           MI(T,R)
%  rho          joint density estimator
%  dD           dpsi*drrho
%  drho         derivative of joint density estimator
%  d2psi        log2(length(Tc))*sqrt(nT*nR);
%
% see also distances/contents
%==============================================================================

function [Dc,rho,dD,drho,d2psi] = MImex(Tc,Rc,omega,m,varargin);

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

% initialize the output variables
Dc  = []; dD = []; rho  = []; drho = []; d2psi = [];

doDerivative = (nargout > 3);  % boolean for derivative computation

%setup default parameter for joint density estimator
tol  = 1e-6;  % tolerance for logarithm in entropy
minT = 0;     % (assumed) smallest value in T
maxT = 256;   % (assumed) largest value in T
nT   = 32;    % number of "bins" for T
minR = 0;     % (assumed) smallest value in R
maxR = 256;   % (assumed) largest value in R
nR   = 30;    % number of "bins" for R

% overwrite default parameter
for k=1:2:length(varargin), % overwrite defaults
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

widthT = (maxT-minT)/nT;
widthR = (maxR-minR)/nR;

% to ensure nothing is missed at the boundary, two artifical bins are added
minT  = minT-2*widthT;  maxT  = maxT+2*widthT;   nT = nT+4;
minR  = minR-2*widthR;  maxR  = maxR+2*widthR;   nR = nR+4;

% compute the Parzen-window estimator for the density and its derivative
% rho  is (nT*nR) x 1
% drho is (nT*nR) x length(Tc), large but sparse
if ~doDerivative,
    try
        rho = rhoSplineC(Tc,Rc,minT,maxT,nT,minR,maxR,nR);
    catch err
        FAIRerror(err);
    end
else
    try
        [rho,drho] = rhoSplineC(Tc,Rc,minT,maxT,nT,minR,maxR,nR);
    catch err
        FAIRerror(err);
    end
    % Attention: This is a hack!
    % Problem: row-indices of sparse matrix have to be ordered ascending
    %          in Matlab.
    %
    % The two lines reorder the sparse matrix drho to make it a valid
    % sparse matrix
    [i,j,s]=find(drho);
    drho = sparse(i,j,s,size(drho,1),size(drho,2));
end;

% reformat rho for efficient computation of marginal densities
% rhoT and rhoR and undo formating
% note, summation can be described via matrix muliplication:
% rhoT = ST*rho, rhoR = SR*rho, where
% ST   = kron( eye(nT,nT), ones(1,nR) ) and
% SR   = kron( ones(1,nT), eye(nR,nR) );
% this is used for the computation of the derivative

rho  = reshape(rho,nT,nR);
rhoT = sum(rho,2);
rhoR = sum(rho,1)';
rho  = rho(:);

% compute MI
Dc = rhoT'*log(rhoT+tol) + rhoR'*log(rhoR+tol) - rho'*log(rho+tol);

if ~doDerivative, return; end;

% build the matrices ST, SR
ST    = sparse(kron(ones(1,nR),speye(nT,nT)));
SR    = sparse(kron(speye(nR,nR),ones(1,nT)));

dpsi  = (log(rhoT+tol)+rhoT./(rhoT+tol))'*ST ...
    +     (log(rhoR+tol)+rhoR./(rhoR+tol))'*SR ...
    -     (log(rho +tol)+rho ./(rho +tol))';

% apply chain rule for MI = psi(rho(T))
dD    = dpsi * drho;

if nargout < 4, return;  end;

% drho  = spdiags(sum(drho,1)',0,length(Tc),length(Tc));
d2psi = log2(length(Tc))*sqrt(nT*nR);
%==============================================================================

