%==============================================================================
% This code is part of the VAMPIRE app for the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR/VAMPIRE.m 
%==============================================================================
% VAMPIRE - Variational Algorithm for Mass-Preserving Image REgistration
% 
% 
%         xxxxxxxx                                           xxxx          
%      xxxxxxxxxxxxxx                                     xxxxxxxxxxx      
%    xxxxxxxxxxxxxxxxx                                 xxxxxxxxxxxxxxxx    
%   xxxxxxxxxxxxxxxxxxxx                             xxxxxxxxxxxxxxxxxxxx  
%  xxxxxxxxxxxxxxxxxxxxxxx                          xxxxxxxxxxxxxxxxxxxxxx 
%       xxxxxxxxxxxxxxxxxxx                       xxxxxxxxxxxxxxxxxxxxx    
%         xxxxxxxxxxxxxxxxxxx     xxx   xxx      xxxxxxxxxxxxxxxxxxx       
%         xxxxxxxxxxxxxxxxxxxx    xxxxxxxxx    xxxxxxxxxxxxxxxxxxxx        
%          xxxxxxxxxxxxxxxxxxxxx  xxxxxxxxx  xxxxxxxxxxxxxxxxxxxxx         
%          xxxx  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx         
%                   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx       xx         
%                     xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx                    
%                      xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx                     
%                      x       xxxxxxxxxxxxxxxx    xxx                     
%                               xxxxxxxxxxxxxx                             
%                               xxx  xxxxxxxx                              
%                               x          xx                              
% 
% 
% function [Jc,para,dJ,H] = VAMPIREPIRobjFctn(T,Rc,omega,m,beta,M,wRef,xc,wc)
% 
% Objective function for parametric mass-preserving transformations
% (see VAMPIRENPIRobjFctn.m for non-parametric MP models)
%
% (1)   J(yc) = J(wc) = D(T(Y(wc)) * det(DY(wc)), R) + alpha * S(wc)
% 
% General outline of the VAMPIRE model:
% -------------------------------------
% The idea behind mass-preserving (MP) registration is that a
% transformation yc applied to an image T must not increase or decrease
% the mass of the image (e.g. in density based modalities like PET).
% With the integration by substitution theorem we get for an image T (y
% diffeomorphic and orientation preserving)
%
% (2)   \int_omega T(x) dx = \int_y(omega) T(y) det(Dy) dx,
% 
% where \int_omega T(x) dx defines the mass (or total amount of intensity)
% of the image T.
%
% Hence, in addition to interpolating template T at the deformed grid y we
% have to apply an intensity modulation given by the Jacobian determinant
% det(Dy) of y.
%
% Example - When for instance a tissue with given amount of signal is
%           stretched such that its volume doubles, the intensity
%           corresponding to the signal should be divided by two.
% 
% 
% Special Requirements:
% ---------------------
% Because of (2) the transformation must be diffeomorphic and orientation
% preserving. Both properties can be achieved by starting with an
% orientation preserving initial guess (e.g. identity) and guaranteeing
% that the volumes of any volume element remains positive at all times.
% To this end (and because of the polyconvex structure of the data term
% [Note D depends on det(Dy)]) we employ the hyper elastic regularization
% functional with projected Armijo line search.
%
%
% Input:
% ------
%   T     - coefficients for interpolation of T
%   Rc    - interpolated R
%   omega - representation of computational domain
%   m     - discretization size
%   beta  - adding beta*I to approximation of Hessian
%   M     - value(s) for regularization of spline coefficients
%           (for M=0, hyper elastic regularization is performed)
%   wRef  - reference parameters
%   xc    - discretization of Omega
%   wc	  - current parameters
%
% Output:
% -------
%   Jc    - current function value J(wc)
%   para  - struct {Tc=T(y(wc)), Rc, omega, m, yc=y(wc,xc), Jc}, for plots
%   dJ    - gradient of J
%   H     - approximation to Hessian of J
%
% Please cite:
%
% @article{GigengackEtAl2012,
%          author = {Gigengack, F and Ruthotto, L and Burger, M and Wolters, C H and Jiang, Xiaoyi and Schafers, K P},
%          title = {{Motion Correction in Dual Gated Cardiac PET Using Mass-Preserving Image Registration}},
%          journal = {Medical Imaging, IEEE Transactions on},
%          year = {2012},
%          volume = {31},
%          number = {3},
%          pages = {698--712},
%          keywords = {Image Registration},
%          doi = {10.1109/TMI.2011.2175402},
%          }

function [Jc,para,dJ,H] = VAMPIREPIRobjFctn(T,Rc,omega,m,beta,M,wRef,xc,wc)

para = []; dJ = []; H = [];
persistent P

% if wc is not an input argument, reports status
if ~exist('wc','var') || isempty(wc),
    if nargout == 1, Jc = 'NPIR';  return; end;
    % report current settings
    dimstr  = @(m) sprintf('[%s]',sprintf(' %d',m));
    wc      = trafo('w0');
    fprintf('Mass-Preserving Parametric Image Registration:');
    fprintf('    J(wc)=D(T(y(wc))*det(Dy(wc)),R) + (wc-wRef)''*M*(wc-wRef) != min\n');
    fprintf('  %20s : %s\n','m',dimstr(m));
    fprintf('  %20s : %s\n','omega',dimstr(omega));
    fprintf('  %20s : %s\n','IMAGE MODEL',imgModel);
    fprintf('  %20s : %s\n','DISTANCE',distance);
    fprintf('  %20s : %s\n','TRAFO',trafo);
    fprintf('  %20s : %s\n','length(wc)',num2str(length(wc)));
    % if ~isempty(M),
    fprintf('  %20s : %s\n','M',sprintf('is %d-by-%d',size(M)));
    % end;
    fprintf('  %20s : %s\n','beta',num2str(beta)); 
    Jc = wc; % return starting guess
    return
end

hd = prod((omega(2:2:end)-omega(1:2:end))./m);
doDerivative = (nargout>2);
matrixFree = regularizer('get','matrixFree');

% transfer to staggered grids and compute cell centered derivatives
P = gridInterpolation(P,omega,m);

% compute transformation, distance, and regularization and combine these
trafo('set', 'm', m+1, 'omega', omega);
[yc,dy] = trafo(wc,xc,'doDerivative',doDerivative);

% intensity modulation
[Jac,dJac] = geometry(yc,m,'Jac','doDerivative',doDerivative,'omega',omega);

% interpolation and intensity modulation
[Tc,dT] = imgModel(reshape(T,m),omega,center(yc,m),'doDerivative',doDerivative,'matrixFree',matrixFree);
Tcmod = Tc .* Jac;

% Note: d1D = dD and d2D = -dD as SSD(I1,I2) ~ (I1-I2)^2/2
[Dc,~,dD,dres,d2psi] = distance(Tcmod,Rc,omega,m,'doDerivative',doDerivative);

% add regularization
if isempty(M) || (all(size(M)==1) && (M==0))
    if strcmp(trafo, 'splineTransformation2D') ...
            || strcmp(trafo, 'splineTransformation2Dsparse') ...
            || strcmp(trafo, 'splineTransformation3Dsparse') ...
            || strcmp(trafo, 'affine3D') ...
            || strcmp(trafo, 'affine3Dsparse')
        [Sc,dS,d2S] = regularizer(yc-getNodalGrid(omega,m),omega,m,'doDerivative',doDerivative);
    else
        d2S = 0;
        dS  = 0;
        Sc  = 0;
        M   = 0;
    end
else
    if isempty(wRef), wRef = trafo('w0'); end
    dS  = (wc-wRef)'*M;
    Sc  = 0.5*dS*(wc-wRef);
end
% joint functional
Jc = Dc + Sc;

% collect variables for plots
para = struct('Tc',Tcmod,'Rc',Rc,'omega',omega,'m',m,'yc',yc,'Jc',Jc);

if ~doDerivative, return; end

if isnumeric(dy) % matrix based code
	%      Tcmod    = T(y) * Jac(y)
	% ==> dTcmoddy  = Jac * dT*P + T(y) * dJac   (product rule)
    dTcmoddy = sdiag(Jac)*dT*P + sdiag(Tc)*dJac;
    
    if isempty(M) || (all(size(M)==1) && (M==0))
        % joint functional
        dJ = (dD*dTcmoddy + dS) * dy;

        if nargout<4, return; end
        % approximation to Hessian
        dr = dres*dTcmoddy;
        H  = dy'*(dr'*d2psi*dr + d2S)*dy;
    else
        % joint functional
        dJ = dD*dTcmoddy * dy + dS;

        if nargout<4, return; end
        % approximation to Hessian
        dr = dres*dTcmoddy*dy;
        H  = dr'*d2psi*dr + M + beta*speye(length(wc));
    end
elseif isstruct(dy) % matrix free code
    if isempty(M) || (all(size(M)==1) && (M==0)) % hyper elastic reg.
        % joint functional
        dJ = dy.Qadjoint(P(vecXmat1(dD,Jac,dT)) ...
            + dJac.dJacadj(yc,m,dD'.*Tc) + dS(:))';
        
        if nargout<4, return; end
        % store information required for approximation of Hessian later
        H.omega         = omega;
        H.m             = m;
        H.d2D.howdTcmod = 'sdiag(Jac)*dT*P + sdiag(Tc)*dJac';
        H.d2D.how       = 'P''* (dres*dTcmod)''*(dres*dTcmod) *P';
        H.d2D.P         = P;
        H.d2D.Tc        = Tc;
        H.d2D.dT        = dT;
        H.d2D.Jac       = Jac;
        H.d2D.dJac      = dJac.dJac;
        H.d2D.dJacadj   = dJac.dJacadj;
        H.d2D.dres      = dres;
        H.d2D.w2y       = @(wc) dy.Q(wc);
        H.d2D.y2w       = @(wc) dy.Qadjoint(wc);
    
        H.d2S         = d2S;
        
    else % spline coefficient regularization
        % joint functional
        dJ = dy.Qadjoint(P(vecXmat1(dD,Jac,dT)) ...
            + dJac.dJacadj(yc,m,dD'.*Tc))' + dS;
        
        if nargout<4, return; end
        dTcmod    = @(x) vecXmat2(Jac,dT,P(x)) + Tc .* dJac.dJac(yc,m,x);
        dTcmodadj = @(x) P(vecXmat2adj(dT,Jac,x)) ...
                       + dJac.dJacadj(yc,m,Tc.*x);
        operator  = @(x) dy.Qadjoint(hd * dTcmodadj(dres' * ...
                        (dres * dTcmod(x))));
        H = @(x) operator(dy.Q(x)) + M*x + beta;
    end
elseif iscell(dy) % cell mode
    if matrixFree % matrix free code
        % joint functional
        if isempty(M) || (all(size(M)==1) && (M==0))
            dJ = kronIdy((P(vecXmat1(dD,Jac,dT)) + dJac.dJacadj(yc,m,dD'.*Tc))' + dS,dy);
        else
            dJ = kronIdy((P(vecXmat1(dD,Jac,dT)) + dJac.dJacadj(yc,m,dD'.*Tc))',dy) + dS;
        end

        if nargout<4, return; end
        dTcmod    = @(x) vecXmat2(Jac,dT,P(x)) + Tc .* dJac.dJac(yc,m,x);
        dTcmodadj = @(x) P(vecXmat2adj(dT,Jac,x)) ...
            + dJac.dJacadj(yc,m,Tc.*x);
        operator  = @(x) kronIdy((hd * dTcmodadj(dres' * ...
            (dres * dTcmod(x))))',dy);
        H = @(x) operator(w2y(dy,x)')' + M*x + beta;
    else % matrix based code
        dTcmod = sdiag(Jac)*dT*P + sdiag(Tc)*dJac;
        % joint functional
        dJ = kronIdy(dD*dTcmod,dy) + dS;

        if nargout<4, return; end
        % approximation to Hessian
        dr = kronIdy(dres*dTcmod,dy);
        H  = dr'*d2psi*dr + M + beta*speye(length(wc));
    end
else
    error('Unknown format of dy!');
end

function P = gridInterpolation(P,omega,m)
switch regularizer
    case 'mbHyperElastic' 
        if size(P,1) ~= length(omega)/2 * prod(m)
            P = nodal2center(m);
        end
    case 'mfHyperElastic'
        P = @(y) nodal2center(y,m);
	otherwise
	    error('VAMPIRE requires Hyperelastic Regularization!');
end

function p = vecXmat1(dD,Jac,dT)
% implementation of vector' * matrix product dD * [Jac] * dT
% where
%
% size(dD)  = [1, nVol * prod(m)]
% size(Jac) = [prod(m),1]
% size(dI)  = [prod(m),dim,nVol]
dim = size(dT,2);
p   = zeros(length(Jac),dim);
for d=1:dim,
    p(:,d) = dD(:).*Jac(:).*dT(:,d);
end
p = reshape(p,1,[]);

function p = vecXmat2(Jac,dI,x)
% We aim to compute diag(Jac) * [dI] * x
%
% size(Jac) = [n, 1]
% size(dI) = [nVol * n , dim * n]
% size(x)  = [dim*n,1]
dim  = size(dI,2);
nVol = size(dI,3);
n    = length(Jac);

x = reshape(x,[],dim);
p = zeros(n,nVol);
for vol=1:nVol,
    for d=1:dim,
        p(:,vol) = p(:,vol) + dI(:,d,vol) .* x(:,d);
    end
    p(:,vol) = p(:,vol).*Jac;
end
p = p(:);

function p = vecXmat2adj(dI,Jac,x)
% We aim to compute dI' * [Jac] * x
%
% size(Jac) = [n,1]
% size(dI)  = [nVol*n, dim*n]
% size(x)   = [nVol*n,1]
dim  = size(dI,2);
nVol = size(dI,3);
n = length(Jac);

x = reshape(x,[],nVol);
p = zeros(n,dim);
for vol=1:nVol,
    for d=1:dim,
        p(:,d) = p(:,d) + dI(:,d,vol) .* x(:,vol);
    end
end
for d=1:dim,
    p(:,d) = p(:,d).*Jac;
end
p = p(:);

function a = sdiag(a)
% shortcut for sparse diagonal matrices
a = spdiags(reshape(a,[],1),0,length(a),length(a));

function dd = kronIdy(Df,Dy)
% It is assumed that Q = I_d \otimes Dy{1},
%
% [df1,...dfd]*| Q       | = [df1*Q,...,dfd*Q]
%              |  \      |
%              |        Q|

m = size(Df,2);
n = size(Dy{1},1);
switch m/n
    case 2, dd = [Df(:,1:n)*Dy{:},...
                  Df(:,n+1:end)*Dy{:}];
    case 3, dd = [Df(:,1:n)*Dy{:},...
                  Df(:,n+1:2*n)*Dy{:},...
                  Df(:,2*n+1:end)*Dy{:}];
    otherwise, error('nyi');
end

function dd = w2y(Q,w)
m = size(Q{:},2);
w = reshape(w,m,[]);
switch size(w,2)
    case 2, dd = [Q{:}*w(:,1);Q{:}*w(:,2)];
    case 3, dd = [Q{:}*w(:,1);Q{:}*w(:,2);Q{:}*w(:,3)];
    otherwise, error('nyi');
end