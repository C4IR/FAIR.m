%==============================================================================
% This code is part of the Finite Element Method app for the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR/FAIRFEM 
%==============================================================================
%
% function [Jc,para,dJ,H] = FEMPIREobjFctn(T,Rc,omega,m,yRef,yc,xc,tri)
%
% Objective Function for Finite Element Mass-Preserving Image Registration
%
% computes J(yc) = SSD(T(yc)*det(Dy), Rc) +alpha*regularizer(yc - xc), where
%
%  Tc,Rc      = T(yc) = imgModel(T,omega,P*yc), P barycentric interpolation of yc
%  SSD(Tc,Rc) = the usual SSD ( mid-point rule on the tetrahedra)
%  S(uc)      = regularizer(uc,omega,m), see elasticFEM and hyperElasticFEM
%
% Modes
%   FEMPIREobjFctn      - displays the help, runs minimal example
%   FEMPIREobjFctn([])  - reports the current setting of the distance and regularization
%   [Jc,para,dJ,H]      = FEMPIREobjFctn(T,Rc,omega,m,yRef,yc,xc,tri)
%                       - evaluates the objective function
%
% References:
%
% [1] Gigengack, F., Ruthotto, L., Burger, M., Wolters, C. H., Jiang, X., & Schafers, K. P. (2012). 
%     Motion Correction in Dual Gated Cardiac PET Using Mass-Preserving Image Registration. 
%     Medical Imaging, IEEE Transactions on, 31(3), 698?712. http://doi.org/10.1109/TMI.2011.2175402
% [2] Ruthotto, L., & Modersitzki, J. (2015). 
%     Non-linear Image Registration. 
%     In Handbook of Mathematical Methods in Imaging (pp. 2005?2051). http://doi.org/10.1007/978-1-4939-0790-8_39
%
% Input
%   T     - data for template image Tc = imgModel(T,omega,yc)
%   Rc    - reference image on barycenters Rc = imgModel(R,omega,P*xc)
%   Mesh  - representation of computational mesh
%   yRef  - coefficients of reference transformation, Sc = regularizer(yc-yRef,omega,m)
%   yc    - current coefficients
%
% Output
%   Jc    - function value
%   para  - parameter for plots
%   dJ    - gradient
%   H     - GaussNewton style approximation to Hessian (either explicit or as operator)
%
% see also VAMPIRENPIRobjFctn (version for nodal grid discretization)
%==============================================================================
function [Jc,para,dJ,H] = FEMPIREobjFctn(T,Rc,Mesh,yRef,yc)

if nargin == 0,
    help(mfilename); runMinimalExample; return;
elseif ~exist('yc','var') || isempty(yc),
    if nargout == 1, Jc = 'FEMPIRE';  return; end;
    dimstr  = @(m) sprintf('[%s]',sprintf(' %d',m));
    fprintf('FEMPIRE - Finite Element Mass-Preserving Image REgistration\n');
    v = @(str) regularizer('get',str); % get regularization configuration like grid
    
    fprintf('  J(yc) = SSD(T(yc)*detDy,R) + alpha*S(yc-yReg) != min\n');
    fprintf('  %20s : %s\n','IMAGE MODEL',imgModel);
    fprintf('  %20s : %s\n','DISTANCE','SSD');
    fprintf('  %20s : %s\n','REGULARIZER',regularizer);
    fprintf('  %20s : %s\n','alpha,mu,lambda',num2str([v('alpha'),v('mu'),v('lambda')]));
    fprintf('  %20s : %s\n','MATRIX FREE',int2str(v('matrixFree')));
    fprintf('  %20s : %s\n','GRID','Finite Element Grid');
    fprintf('  %20s : %s\n','TETRAHDRA',num2str(size(Mesh.tri,1)));
    fprintf('  %20s : %s\n','NODES',num2str(size(Mesh.xn,1)));
    fprintf('  %20s : %s\n','m',dimstr(Mesh.m));
    fprintf('  %20s : %s\n','omega',dimstr(Mesh.omega));
    return;
end;
dim   = Mesh.dim;
omega = Mesh.omega;
m     = Mesh.m;

yc = reshape(yc,[],dim);
yRef = reshape(yRef,[],dim);
doDerivative = (nargout>2);            % flag for necessity of derivatives
matrixFree = regularizer('get','matrixFree');
% do the work ------------------------------------------------------------

% define barycentric interpolation

% compute weights for mid-point quadrature
vol = Mesh.vol;

% compute interpolated image and derivative
[Tc,dT] = imgModel(T,omega,Mesh.mfPi(yc,'C'),'doDerivative',doDerivative,'matrixFree',matrixFree);
% apply intensity modulation
[det,dDet]  = getDeterminant(Mesh,yc,'doDerivative',doDerivative);
Tcmod   = Tc .* det;

% compute SSD distance
rc = Tcmod-Rc;                     % the residual
Dc = 0.5*(rc'*sdiag(vol)*rc);       	    % the SSD
dres = 1;
dD =  rc'*sdiag(vol)*dres;
d2psi = sdiag(vol);

% compute regularizer, note that xc,yRef and triangulation must be supplied
[Sc,dS,d2S] = regularizer(yc(:)-yRef(:),yRef(:),Mesh,'doDerivative',doDerivative);

% evaluate joint function
Jc = Dc + Sc;

% store intermediates for outside visualization
para = struct('Tc',Tcmod,'Rc',Rc,'m',m,'yc',yc,'Mesh',Mesh,'Jc',Jc,'Dc',Dc,'Sc',Sc);

if ~doDerivative, return; end;

if not(regularizer('get','matrixFree')),
    % matrix based mode, business as usual
    %
    %      Tcmod(y)  = T(y) * Jac(y)
    % ==> dTcmod(y)  = Jac * dT * P + T(y) * dJac   (product rule)
    P = kron(speye(dim),Mesh.PC);
    dTcmod = sdiag(det)*dT*P + sdiag(Tc)*dDet;
    % chain rule: D(y) = D(I1,I2) = D(Tcmod(y),R) ==> dD = d_I1 D * dTcmod
    dr = dres*dTcmod;
    dD = dD*dTcmod;
    
    dJ = dD + dS;
    H  = dr'*d2psi*dr + d2S;
else
    % derivatives rather explicit
    %P  = @(x) reshape(Mesh.mfPi(reshape(x,[],dim),'C'),[],1);
    %dTcmod = P((sdiag(det)*dT)')' + sdiag(Tc)*dDet;    
    p = (dD(:).*det).*dT;
    dD = reshape(Mesh.mfPi(p,'C'),1,[]) + (dDet.dDetadj(dD'.*Tc))';
    
    dJ = dD + dS;
    
    % approximation to d2D in matrix free mode
    % d2D   = P'*dr'*d2psi*dr*P
    % P and P' are operators matrix free
    H.Mesh      = Mesh;
    H.omega     = omega;
    H.m         = m;
    H.d2D.how   = 'P''*dr''*d2psi*dr*P';
    H.d2D.P     = @(x) reshape(Mesh.mfPi(x,'C'),[],1);
    H.d2D.Tc    = Tc;
    H.d2D.dT    = dT;
    H.d2D.Jac   = det;
    H.d2D.dJac  = dDet.dDet;
    H.d2D.dJacadj = dDet.dDetadj;
    H.d2D.diag  = @(Mesh,yc) getDiag(Mesh,det,Tc,dT,yc);
    H.d2D.dres  = dres;
    H.d2D.d2psi = d2psi;
    H.solver    = d2S.solver;
    
    H.d2S = d2S;
end;


% shortcut for sparse diagonal matrix
function A = sdiag(v)
A = spdiags(v(:),0,numel(v),numel(v));

% diagonal of hessian matrix
% sum(sdiag(vol)dTcmod.^2,1)
% where dTcmod = sdiag(det)*dT*P + sdiag(Tc)*dDet;
% implemention:  sum(sdiag(vol)*(sdiag(det)*dT*P(.^2 +
% sdiag(vol)*(sdiag(Tc)*dDet).^2,1)
function D = getDiag(Mesh,det,Tc,dT,yc)

vol = Mesh.vol;
dim = Mesh.dim;

if dim == 2

    % first term
    D1 = Mesh.mfPi(vol.*(det.*dT).^2,'C');
    D1 = D1(:)./3;

    % second term
    dphi = Mesh.dphi;    
    dphi1 = dphi{1}; dphi2 = dphi{2}; dphi3 = dphi{3};
    
    % compute diagonal of volume regularizer
    % MB : dDet = [sdiag(By(:,4)), sdiag(-By(:,3)) , sdiag(-By(:,2)), sdiag(By(:,1))] * B;
    [dx1, dx2] = getGradientMatrixFEM(Mesh,1);
    yc = reshape(yc,[],2);
    By = [dx1.D(yc(:,1))  dx2.D(yc(:,1)) dx1.D(yc(:,2)) dx2.D(yc(:,2))];

    % get boundaries
    Dxi = @(i,j) Mesh.mfPi(vol.*(Tc.*By(:,j).*dphi1(:,i)).^2,1) ... 
        + Mesh.mfPi(vol.*(Tc.*By(:,j).*dphi2(:,i)).^2,2) ...
        + Mesh.mfPi(vol.*(Tc.*By(:,j).*dphi3(:,i)).^2,3); % diagonal of Dxi'*Dxi 

    Dxy = @(i,j,k,l) Mesh.mfPi(vol.*(Tc.*By(:,j).*dphi1(:,i)).*(Tc.*By(:,l).*dphi1(:,k)),1) ... % byproduct terms for verifications
        + Mesh.mfPi(vol.*(Tc.*By(:,j).*dphi2(:,i)).*(Tc.*By(:,l).*dphi2(:,k)),2) ...
        + Mesh.mfPi(vol.*(Tc.*By(:,j).*dphi3(:,i)).*(Tc.*By(:,l).*dphi3(:,k)),3); % diagonal of Dxi'*Dxi 

    D2 = [Dxi(1,4)+ Dxi(2,3) - 2*Dxy(1,4,2,3) ;Dxi(1,2)+ Dxi(2,1) - 2*Dxy(1,2,2,1)];

    % third byproduct term
    D3 = [(Mesh.mfPi(vol.*((1/3)*det.*dT(:,1)).*Tc.*(By(:,4).*dphi1(:,1) - By(:,3).*dphi1(:,2)),1)) + ...
          (Mesh.mfPi(vol.*((1/3)*det.*dT(:,1)).*Tc.*(By(:,4).*dphi2(:,1) - By(:,3).*dphi2(:,2)),2)) + ...
          (Mesh.mfPi(vol.*((1/3)*det.*dT(:,1)).*Tc.*(By(:,4).*dphi3(:,1) - By(:,3).*dphi3(:,2)),3)) ;...

        (Mesh.mfPi(vol.*((1/3)*det.*dT(:,2)).*Tc.*(By(:,1).*dphi1(:,2) - By(:,2).*dphi1(:,1)),1)) + ...
        (Mesh.mfPi(vol.*((1/3)*det.*dT(:,2)).*Tc.*(By(:,1).*dphi2(:,2) - By(:,2).*dphi2(:,1)),2)) + ...
        (Mesh.mfPi(vol.*((1/3)*det.*dT(:,2)).*Tc.*(By(:,1).*dphi3(:,2) - By(:,2).*dphi3(:,1)),3))];

    D = D1 + D2 + 2*D3;
    
else
    
    % first term
    D1 = Mesh.mfPi(vol.*(det.*dT).^2,'C');
    D1 = D1(:)./4;
    
    dphi = Mesh.dphi;    
    dphi1 = dphi{1}; dphi2 = dphi{2}; dphi3 = dphi{3}; dphi4 = dphi{4};
    
    [cof, dCof] = cofactor3D([],Mesh,0,1);
    yc = reshape(yc,[],3);
    dx1 = Mesh.mfdx1; dx2 = Mesh.mfdx2; dx3 = Mesh.mfdx3;
    %det = dx1.D(yc(:,1)).*cof{1}(yc) + dx2.D(yc(:,1)).*cof{2}(yc) + dx3.D(yc(:,1)).*cof{3}(yc);
    
    % get boundaries
    Dxi = @(i,j) Mesh.mfPi(vol.*(Tc.*cof{j}(yc).*dphi1(:,i)).^2,1) ... 
        + Mesh.mfPi(vol.*(Tc.*cof{j}(yc).*dphi2(:,i)).^2,2) ...
        + Mesh.mfPi(vol.*(Tc.*cof{j}(yc).*dphi3(:,i)).^2,3) ...
        + Mesh.mfPi(vol.*(Tc.*cof{j}(yc).*dphi4(:,i)).^2,4); % diagonal of Dxi'*Dxi 
    
    Dxy = @(i,j,k,l) Mesh.mfPi(vol.*(Tc.*cof{j}(yc).*dphi1(:,i)).*(Tc.*cof{l}(yc).*dphi1(:,k)),1) ... % byproduct terms for verifications
        + Mesh.mfPi(vol.*(Tc.*cof{j}(yc).*dphi2(:,i)).*(Tc.*cof{l}(yc).*dphi2(:,k)),2) ...
        + Mesh.mfPi(vol.*(Tc.*cof{j}(yc).*dphi3(:,i)).*(Tc.*cof{l}(yc).*dphi3(:,k)),3) ...
        + Mesh.mfPi(vol.*(Tc.*cof{j}(yc).*dphi4(:,i)).*(Tc.*cof{l}(yc).*dphi4(:,k)),4); % diagonal of Dxi'*Dxi 
    
    D21 = Dxi(1,1)+ Dxi(2,2) + Dxi(3,3) + 2*Dxy(1,1,2,2) + 2*Dxy(1,1,3,3) + 2*Dxy(2,2,3,3);
    D22 = Dxi(1,4)+ Dxi(2,5) + Dxi(3,6) + 2*Dxy(1,4,2,5) + 2*Dxy(1,4,3,6) + 2*Dxy(2,5,3,6);
    D23 = Dxi(1,7)+ Dxi(2,8) + Dxi(3,9) + 2*Dxy(1,7,2,8) + 2*Dxy(1,7,3,9) + 2*Dxy(2,8,3,9);
    
    D2 = [D21;D22;D23];
    
    % third term
    D3 = [(Mesh.mfPi(vol.*((1/4)*det.*dT(:,1)).*Tc.*(cof{1}(yc).*dphi1(:,1) + cof{2}(yc).*dphi1(:,2) + cof{3}(yc).*dphi1(:,3)),1)) + ...
          (Mesh.mfPi(vol.*((1/4)*det.*dT(:,1)).*Tc.*(cof{1}(yc).*dphi2(:,1) + cof{2}(yc).*dphi2(:,2) + cof{3}(yc).*dphi2(:,3)),2)) + ...
          (Mesh.mfPi(vol.*((1/4)*det.*dT(:,1)).*Tc.*(cof{1}(yc).*dphi3(:,1) + cof{2}(yc).*dphi3(:,2) + cof{3}(yc).*dphi3(:,3)),3)) + ...
          (Mesh.mfPi(vol.*((1/4)*det.*dT(:,1)).*Tc.*(cof{1}(yc).*dphi4(:,1) + cof{2}(yc).*dphi4(:,2) + cof{3}(yc).*dphi4(:,3)),4));...

        (Mesh.mfPi(vol.*((1/4)*det.*dT(:,2)).*Tc.*(cof{4}(yc).*dphi1(:,1) + cof{5}(yc).*dphi1(:,2) + cof{6}(yc).*dphi1(:,3)),1)) + ...
        (Mesh.mfPi(vol.*((1/4)*det.*dT(:,2)).*Tc.*(cof{4}(yc).*dphi2(:,1) + cof{5}(yc).*dphi2(:,2) + cof{6}(yc).*dphi2(:,3)),2)) + ...
        (Mesh.mfPi(vol.*((1/4)*det.*dT(:,2)).*Tc.*(cof{4}(yc).*dphi3(:,1) + cof{5}(yc).*dphi3(:,2) + cof{6}(yc).*dphi3(:,3)),3)) + ...
        (Mesh.mfPi(vol.*((1/4)*det.*dT(:,2)).*Tc.*(cof{4}(yc).*dphi4(:,1) + cof{5}(yc).*dphi4(:,2) + cof{6}(yc).*dphi4(:,3)),4));...
        
        (Mesh.mfPi(vol.*((1/4)*det.*dT(:,3)).*Tc.*(cof{7}(yc).*dphi1(:,1) + cof{8}(yc).*dphi1(:,2) + cof{9}(yc).*dphi1(:,3)),1)) + ...
        (Mesh.mfPi(vol.*((1/4)*det.*dT(:,3)).*Tc.*(cof{7}(yc).*dphi2(:,1) + cof{8}(yc).*dphi2(:,2) + cof{9}(yc).*dphi2(:,3)),2)) + ...
        (Mesh.mfPi(vol.*((1/4)*det.*dT(:,3)).*Tc.*(cof{7}(yc).*dphi3(:,1) + cof{8}(yc).*dphi3(:,2) + cof{9}(yc).*dphi3(:,3)),3)) + ...
        (Mesh.mfPi(vol.*((1/4)*det.*dT(:,3)).*Tc.*(cof{7}(yc).*dphi4(:,1) + cof{8}(yc).*dphi4(:,2) + cof{9}(yc).*dphi4(:,3)),4))];

    D = D1 + D2 + 2*D3;
end
    

function runMinimalExample
setup2DGaussianData

imgModel('reset','imgModel','splineInter','regularizer','moments','theta',1e0);
regularizer('reset','regularizer','mbHyperElasticFEM','alpha',1e1,'mu',1,'lambda',0);
level = 4; omega = ML{level}.omega; m = ML{level}.m;
[T,R] = imgModel('coefficients',ML{level}.T,ML{level}.R,omega(1,:),'out',0);

Mesh  = TriMesh1(omega,m);

% interpolate reference
Rc = imgModel(R,omega,Mesh.mfPi(Mesh.xn,'C'));

yc = Mesh.xn + 1e-2 * randn(size(Mesh.xn));
fctn = @(yc) feval(mfilename,T,Rc,Mesh,Mesh.xn(:),yc(:));

[J,para,dJ,H] = fctn(yc);
checkDerivative(fctn,yc(:));
