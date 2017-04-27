%==============================================================================
% This code is part of the Finite Element Method app for the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR/FAIRFEM 
%==============================================================================
%
% function [Sc,dS,d2S] = elasticFEM(uc,xc,tri,varargin)
%
% linear elastic regularizer for tetrahedral finite element discretization
% with linear basis functions
%
% S(u) = \int || B * u  ||^2 dx
%
% where B is the navier lame operator
%
% Input:
%   uc    - coefficients of current displacement
%   yRef  - coefficients of reference transformation, yc = uc + yRef
%   Mesh  - representation of mesh
%
% Output:
%    Sc    - hyperelastic energy
%    dS    - first derivative
%    d2S   - Hessian
%
% see also elastic.m (staggered grid version)
%==============================================================================

function [Sc,dS,d2S] = elasticFEM(uc,~,Mesh,varargin)

if nargin == 0, help(mfilename); runMinimalExample;return;end

persistent A alphaOld
if ~exist('A','var'),        A        = []; end;
if ~exist('alphaOld','var'), alphaOld = []; end;

% default parameter
matrixFree  = 0;
alpha       = 1;
mu          = 1;
lambda      = 0;

for k=1:2:length(varargin), % overwrites default parameter
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

if not(matrixFree), % matrix-based
    build = isempty(alphaOld) || any(alphaOld~=[alpha mu lambda]) ||...
        isempty(A) || size(A,2)~=numel(uc);
    if build,
        alphaOld = [alpha mu lambda];
        A = getElasticMatrixFEM(Mesh,mu,lambda);
        A = alpha * (A'* A);
    end;
    dS  = uc'*A;
    Sc  = 0.5*dS*uc;
    d2S = A;
else % matrix-free
    d2S.regularizer = regularizer;
    d2S.alpha  = alpha;
    d2S.By     = @(u,Mesh) elasticOperator(u,Mesh,mu,lambda,'By');
    d2S.BTy    = @(u,Mesh) elasticOperator(u,Mesh,mu,lambda,'BTy');
    d2S.B      = @(Mesh) getElasticMatrixFEM(Mesh,mu,lambda);
    d2S.diag   = @(Mesh) getDiag(Mesh,mu,lambda,alpha);
    d2S.solver = @FEMMultiGridsolveGN;
    d2S.d2S  = @(uc,Mesh) ...
        alpha * ...
        elasticOperator(elasticOperator(uc,Mesh,mu,lambda,'By'),Mesh,mu,lambda,'BTy');
    
    dS   = d2S.d2S(uc,Mesh)';
    Sc   = .5*dS*uc;
end

function By = elasticOperator(u,Mesh,mu,lambda,flag)
vol  = sqrt(Mesh.vol);
a    = sdiag(sqrt(mu)*vol);
b    = sdiag(sqrt(mu+lambda)*vol);
flag = [num2str(Mesh.dim) 'D-' flag];
switch flag
    case '2D-By'
        u  = reshape(u,[],2);
        dx1 = Mesh.mfdx1.D ; dx2 = Mesh.mfdx2.D;
        
        By = [...
            a*dx1(u(:,1));...
            a*dx2(u(:,1));...
            a*dx1(u(:,2));...
            a*dx2(u(:,2));...
            b*dx1(u(:,1))+b*dx2(u(:,2))...
            ];
    case '2D-BTy'
        u = reshape(u,[],5);
        dx1 = Mesh.mfdx1.Dadj ; dx2 = Mesh.mfdx2.Dadj;
        
        By = [dx1(a*u(:,1))+dx2(a*u(:,2))+dx1(b*u(:,5));...
            dx1(a*u(:,3))+dx2(a*u(:,4))+dx2(b*u(:,5))];
    case '3D-By'
        u  = reshape(u,[],3);
        dx1 = Mesh.mfdx1.D ; dx2 = Mesh.mfdx2.D; dx3 = Mesh.mfdx3.D;
        
        By = [...
                a*dx1(u(:,1));...
                a*dx2(u(:,1));...
                a*dx3(u(:,1));...
                a*dx1(u(:,2));...
                a*dx2(u(:,2));...
                a*dx3(u(:,2));...
                a*dx1(u(:,3));...
                a*dx2(u(:,3));...
                a*dx3(u(:,3));...
                b*dx1(u(:,1))+b*dx2(u(:,2))+b*dx3(u(:,3))...
            ];
    case '3D-BTy'
        u = reshape(u,[],10);
        dx1 = Mesh.mfdx1.Dadj ; dx2 = Mesh.mfdx2.Dadj; dx3 = Mesh.mfdx3.Dadj;
       
        By = [dx1(a*u(:,1))+dx2(a*u(:,2))+dx3(a*u(:,3))+dx1(b*u(:,10));...
            dx1(a*u(:,4))+dx2(a*u(:,5))+dx3(a*u(:,6))+dx2(b*u(:,10));...
            dx1(a*u(:,7))+dx2(a*u(:,8))+dx3(a*u(:,9))+dx3(b*u(:,10))...
            ];
        
    otherwise
        error('unknown flag: %s',flag);
end


% get diagonal of d2S (interesting in matrix free mode)
function D = getDiag(Mesh,mu,lambda,alpha)
dim = Mesh.dim;
vol = Mesh.vol;
if dim==2
    xn = Mesh.xn;
    x1 =  Mesh.mfPi(xn,1);
    x2 =  Mesh.mfPi(xn,2);
    x3 =  Mesh.mfPi(xn,3);
    e1 =  x1 - x3;
    e2 =  x2 - x3;
    
    % compute gradients of basis functions
    dphi1 =   [ e2(:,2) -e2(:,1)]./[2*vol,2*vol];
    dphi2 =   [-e1(:,2)  e1(:,1)]./[2*vol,2*vol];
    dphi3 = -dphi1 - dphi2;
    
    
    % get boundaries
    Dxi = @(i)   Mesh.mfPi(vol.*dphi1(:,i).^2,1) ...
        + Mesh.mfPi(vol.*dphi2(:,i).^2,2) ...
        + Mesh.mfPi(vol.*dphi3(:,i).^2,3); % diagonal of Dx1'*Dx1
    
    
    
    D1 = (2*mu+lambda)*Dxi(1)+            mu*Dxi(2);
    D2 = mu*Dxi(1)           + (2*mu+lambda)*Dxi(2);
    D = [D1;D2];
else
    xn = Mesh.xn;
    % compute edges
    v1   = Mesh.mfPi(xn,1);
    v2   = Mesh.mfPi(xn,2);
    v3   = Mesh.mfPi(xn,3);
    v4   = Mesh.mfPi(xn,4);
    e1   =  v1 - v4;
    e2   =  v2 - v4;
    e3   =  v3 - v4;
    % compute inverse transformation to reference element
    Ainv =  [
        e2(:,2).*e3(:,3)-e2(:,3).*e3(:,2), ...
        -(e1(:,2).*e3(:,3)-e1(:,3).*e3(:,2)), ...
        e1(:,2).*e2(:,3)-e1(:,3).*e2(:,2), ...
        -(e2(:,1).*e3(:,3)-e2(:,3).*e3(:,1)), ...
        e1(:,1).*e3(:,3)-e1(:,3).*e3(:,1), ...
        -(e1(:,1).*e2(:,3)-e1(:,3).*e2(:,1)), ...
        e2(:,1).*e3(:,2)-e2(:,2).*e3(:,1), ...
        -(e1(:,1).*e3(:,2)-e1(:,2).*e3(:,1)), ...
        e1(:,1).*e2(:,2)-e1(:,2).*e2(:,1), ...
        ];
    detA = e1(:,1).*Ainv(:,1) + e2(:,1).*Ainv(:,2)+ e3(:,1).*Ainv(:,3);
    
    
    % compute gradients of basis functions
    dphi1 =   Ainv(:,1:3)./repmat(detA,[1 3]);
    dphi2 =   Ainv(:,4:6)./repmat(detA,[1 3]);
    dphi3 =   Ainv(:,7:9)./repmat(detA,[1 3]);
    dphi4 =   1 - dphi1 - dphi2 - dphi3;
    
    Dxi = @(i) Mesh.mfPi(vol.*dphi1(:,i).^2,1) +  Mesh.mfPi(vol.*dphi2(:,i).^2,2)  ...
        + Mesh.mfPi(vol.*dphi3(:,i).^2,3) +  Mesh.mfPi(vol.*dphi4(:,i).^2,4);
    
    D1 = (2*mu+lambda)*Dxi(1)+            mu*Dxi(2) + mu *Dxi(3);
    D2 = (2*mu+lambda)*Dxi(2)+            mu*Dxi(1) + mu *Dxi(3);
    D3 = (2*mu+lambda)*Dxi(3)+            mu*Dxi(1) + mu *Dxi(2);
    D = [D1;D2;D3];
end
D = alpha*reshape(D,[],1);


function B = getElasticMatrixFEM(Mesh,mu,lambda)
vol = sqrt(Mesh.vol);
a   = sdiag(sqrt(mu)*vol);
b   = sdiag(sqrt(mu+lambda)*vol);
switch Mesh.dim
    case 2
        dx1 = Mesh.dx1; dx2 = Mesh.dx2;
        B = [a*dx1;a*dx2];
        B = [ blkdiag(B,B); b*dx1, b*dx2];
    case 3
        dx1 = Mesh.dx1; dx2 = Mesh.dx2; dx3 = Mesh.dx3;
        B = [a*dx1;a*dx2;a*dx3];
        B = [ blkdiag(B,B,B); b*dx1, b*dx2, b*dx3];
end

% shortcut for sparse diagonal matrix
function A = sdiag(v)
A = spdiags(v(:),0,numel(v),numel(v));

function runMinimalExample
% ====== 2D
omega = [0,10,0,8]; m = [17,16]; p = [5,6];
w = zeros([p,2]);  w(3,3,1) = 0.06; w(3,4,2) = -0.05;
xn = reshape(getNodalGrid(omega,m),[],2);
yn = splineTransformation2D(w(:),xn(:),'omega',omega,'m',m+1,'p',p,'Q',[]);

Mesh = TriMesh1(omega,m);

regOptn = {'alpha',1,'mu',1,'lambda',0};

fctn = @(y) feval(mfilename,y-xn(:),xn(:),Mesh,regOptn{:},'matrixFree',0);
[Sc, dS, d2S] = fctn(yn(:));
checkDerivative(fctn,yn(:));

[Scmf, dSmf, d2Smf] = feval(mfilename,yn(:)-xn(:),xn(:),Mesh,regOptn{:},'matrixFree',1);

err1 = abs(Sc-Scmf)/abs(Sc);
err2 = norm(dS-dSmf)/norm(dS);

y = randn(5*size(Mesh.tri,1),1); x = rand(numel(xn),1);
err3 = norm(d2S*x-d2Smf.d2S(x,Mesh))/norm(d2S*x);
err4 = abs(y'*d2Smf.By(x,Mesh)-x'*d2Smf.BTy(y,Mesh));
fprintf('%s, MFerror(S,dS,d2S): [%1.2e,%1.2e,%1.2e] OK? %d\n\t\tAdjointError: %1.2e, OK? %d\n',...
    mfilename,err1,err2,err3,max([err1;err2;err3])<1e-14,err4,err4<1e-13);
err5 = norm(full(diag(d2S))-d2Smf.diag(Mesh))/norm(full(diag(d2S)));
fprintf('\t\tDiagErrorMF: %1.2e OK? %d\n',err5,err5<1e-13);

% ====== 3D
omega = [0,10,0,8,0 3]; m = [3,5,6];
Mesh = TetraMesh1(omega,m);

xn = Mesh.xn;
yn = xn + randn(size(xn));

regOptn = {'alpha',1,'mu',1,'lambda',0};

fctn = @(y) feval(mfilename,y-xn(:),xn(:),Mesh,regOptn{:},'matrixFree',0);
[Sc, dS, d2S] = fctn(yn(:));
checkDerivative(fctn,yn(:));

[Scmf, dSmf, d2Smf] = feval(mfilename,yn(:)-xn(:),xn(:),Mesh,regOptn{:},'matrixFree',1);

err1 = abs(Sc-Scmf)/abs(Sc);
err2 = norm(dS-dSmf)/norm(dS);

y = randn(10*size(Mesh.tri,1),1); x = rand(numel(xn),1);
err3 = norm(d2S*x-d2Smf.d2S(x,Mesh))/norm(d2S*x);
err4 = abs(y'*d2Smf.By(x,Mesh)-x'*d2Smf.BTy(y,Mesh));
fprintf('%s, MFerror(S,dS,d2S): [%1.2e,%1.2e,%1.2e] OK? %d\n\t\tAdjointError: %1.2e, OK? %d\n',...
    mfilename,err1,err2,err3,max([err1;err2;err3])<1e-14,err4,err4<1e-13);
err5 = norm(full(diag(d2S))-d2Smf.diag(Mesh))/norm(full(diag(d2S)));
fprintf('\t\tDiagErrorMF: %1.2e OK? %d\n',err5,err5<1e-13);