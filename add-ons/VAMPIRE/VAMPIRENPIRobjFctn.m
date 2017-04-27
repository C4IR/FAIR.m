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
% function [Jc,para,dJ,H] = VAMPIRENPIRobjFctn(T,Rc,omega,m,~,yc)
% 
% Objective function for non-parametric mass-preserving transformations
% (see VAMPIREPIRobjFctn.m for parametric MP models)
% 
% (1)   J(yc) = D(T(yc) * det(Dy), R) + alpha * hyperElastic(yc - xc)
% 
% General outline of the VAMPIRE model:
% -------------------------------------
% The idea behind mass-preserving (MP) registration is that a
% transformation yc applied to an image T must not increase or decrease
% the mass of the image (e.g. in density based modalities like PET).
% With the integration by substitution theorem we get for an image T (y
% diffeomorphic and orientation preserving)
%
% (2)   \int_omega T(x) dx = \int_y(omega) T(y) * det(Dy) dx,
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
%   yc    - current grid
% 
% Output:
% ------
%   Jc    - objective function value
%   para  - struct(Tc,Rc,omega,m,P*yc,Jc), used for plots
%   dJ    - derivative of J, 1-by-length(Y)
%   H     - approximation to Hessian (matrix based or matrix free)
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

function [Jc,para,dJ,H] = VAMPIRENPIRobjFctn(T,Rc,omega,m,~,yc)

Jc = []; para = []; dJ = []; H = [];
persistent P

% if d is not an input argument, reports status
if ~exist('yc','var') || isempty(yc),
    if nargout == 1, Jc = 'NPIR';  return; end;
    dimstr = @(m) sprintf('[%s]',[sprintf('%d',m(1)) sprintf(' %d',m(2:end))]);
    fprintf('VAMPIRE - Variational Algorithm for Mass-Preserving Image REgistration\n');
    v = @(str) regularizer('get',str); % get regularization configuration like grid
    fprintf('  J(yc) = D(T(yc)* det(Dy),R) + alpha * hyperElastic(yc-xc)\n');
    fprintf('  %20s : %s\n','INTERPOLATION',imgModel);
    fprintf('  %20s : %s\n','DISTANCE',distance);
    fprintf('  %20s : %s\n','REGULARIZER',regularizer);
    fprintf('  %20s : [%s]\n','alpha',num2str(v('alpha')));
    fprintf('  %20s : [%s]\n','alphaLength',num2str(v('alphaLength')));
    fprintf('  %20s : [%s]\n','alphaArea',num2str(v('alphaArea')));
    fprintf('  %20s : [%s]\n','alphaVolume',num2str(v('alphaVolume')));
    fprintf('  %20s : %s\n','GRID',v('grid'));
    fprintf('  %20s : %s\n','MATRIX FREE',int2str(v('matrixFree')));
    fprintf('  %20s : %s\n','m',dimstr(m));
    fprintf('  %20s : %s\n','omega',dimstr(omega));
    return
end

hd  = prod((omega(2:2:end)-omega(1:2:end))./m);
doDerivative = (nargout>2);
matrixFree = regularizer('get','matrixFree');

% operator for nodal-->cell centered projections
P = gridInterpolation(P,omega,m);

% compute jacobian determinant for intensity modulation
[Jac,dJac] = geometry(yc,m,'Jac','doDerivative',doDerivative,'omega',omega);

% interpolation and apply intensity modulation
[Tc,dT] = imgModel(reshape(T,m),omega,center(yc,m),'doDerivative',doDerivative,'matrixFree',matrixFree);
Tcmod = Tc .* Jac;

% distance functional
[Dc,~,dD,dres] = distance(Tcmod,Rc,omega,m,'doDerivative',doDerivative);

% regularization functional
[Sc,dS,d2S] = regularizer(yc-getNodalGrid(omega,m),omega,m,'doDerivative',doDerivative);

% joint functional
Jc =  Dc + Sc;

% store intermediates for outside visualization
para = struct('Tc',Tcmod,'Rc',Rc,'omega',omega,'m',m,'yc',center(yc,m),'Jc',Jc);

if ~doDerivative, return; end

if not(matrixFree) % matrix based code
    %      Tcmod(y)  = T(y) * Jac(y)
    % ==> dTcmod(y)  = Jac * dT * P + T(y) * dJac   (product rule)
    dTcmod = sdiag(Jac)*dT*P + sdiag(Tc)*dJac; 
    % chain rule: D(y) = D(I1,I2) = D(Tcmod(y),R) ==> dD = d_I1 D * dTcmod
    dD = dD*dTcmod;
    
    % joint functional
    dJ = dD + dS;
    
    if nargout<4, return; end
    % approximation to Hessian in Gauss-Newton style
    dr = dres*dTcmod;
    H  = hd*(dr'*dr) + d2S;
else % matrix free code
  	% note out = vector'*operator, hence out = (operatorAdjoint(vector))'
    % Hence compute directly dD*dTcmod
    dD = (P(vecXmat(dD,Jac,dT)) + dJac.dJacadj(yc,m,dD'.*Tc))';
    
	% joint functional
    dJ = dD + dS;
    
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
	H.d2D.w2y       = @(x) x;
	H.d2D.y2w       = @(x) x;
 
	H.d2S           = d2S;
    
end

function p = vecXmat(dD,Jac,dT)
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

function a = sdiag(a)
% shortcut for sparse diagonal matrices
a = spdiags(reshape(a,[],1),0,length(a),length(a));

function P = gridInterpolation(P,omega,m)
switch regularizer
    case 'mbHyperElastic',
        if size(P,1) ~= length(omega)/2 * prod(m)
            P = nodal2center(m);
        end
    case {'mfHyperElastic'},
        P = @(y) nodal2center(y,m);
    otherwise
        error('VAMPIRE requires Hyperelastic Regularization!');
end