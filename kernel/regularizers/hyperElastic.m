%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function [S, dS, d2S] = hyperElastic(uc,omega, m, varargin)
%
% computes hyperelastic regularization energy for image registration
%
% Input:
% ------
%    uc         deformation (nodal)
%    omega         spatial domain
%    m             number of discretization points
%    varargin   optional parameter
%
% Output:
% ------
%    S           hyperelastic energy
%    dS          first derivative
%    d2S          approximated Hessian
%  if ~matrixFree,  d2S is sparse matrix; else, d2S is struct endif        
%
% See also
%  @article{2011-BMR,
%    Author = {Burger M., Modersitzki J., Ruthotto L. },
%    Publisher = {University of Muenster},
%    Title = {A hyperelastic regularization energy for image registration},
%    Year = {2011}
%  }
% 
% see also E0_2DCircleToC_hyperElastic
% =============================================================================

function [S, dS, d2S] = hyperElastic(uc, omega, m, varargin)

if nargin == 0, runMinimalExample; return; end;

S = 0; dS = []; d2S = [];
doDerivative = (nargout >1);
% set defaults
alpha       = regularizer('get','alpha');
if isempty(alpha), alpha=1; end 
alphaLength = regularizer('get','alphaLength'); 
if isempty(alphaLength), alphaLength=1; end 
alphaArea   = regularizer('get','alphaArea'); 
if isempty(alphaArea), alphaArea=1;     end 
alphaVolume = regularizer('get','alphaVolume'); 
if isempty(alphaVolume), alphaVolume=1; end;

phiArea = 'phiAreaDoubleWell';

matrixFree = false;
wc = [];
for k=1:2:length(varargin) % overwrites defaults
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end


alphaLength = alpha * alphaLength; 
alphaArea   = alpha * alphaArea; 
alphaVolume = alpha * alphaVolume;

dim = length(omega)/2;
h   = (omega(2:2:end) - omega(1:2:end)) ./ m;
hd  = prod(h);

% note, the input uc is a displacement, while yc=xc+uc is the transformation
xc = getNodalGrid(omega,m);
if not(isempty(wc)),
    xc = trafo(wc,xc);
end
yc       = xc + uc;

% compute reference areas and volumes of standard voxel [0,h1]x[0,h2]x[0,h3]
% note, reference grid is nodal 
mOne  = ones(size(m)); 
voxel = getNodalGrid(reshape([0;1]*h,1,[]),mOne);

%
% first term: define Slength, prepare dSlength and d2Slength
%
% Slength(uc) = .5 * \int norm(gradient uc) dx
if alphaLength~=0
  if not(matrixFree),
    B         = getGradientNodal(omega,m);
    d2Slength = alphaLength * hd * (B' * B);
    dSlength  = d2Slength * uc;
    Slength   = 0.5 * uc' * dSlength ;
  else
    d2Slength = @(x,omega,m) (alphaLength*hd)*...
      lengthOperator(lengthOperator(x,omega,m,'By'),omega,m,'BTy');
    dSlength = d2Slength(uc,omega,m)';
    Slength  = .5*dSlength*uc;
  end
else
  Slength  = 0;
  dSlength = zeros(dim*prod(m+1),1);
  if not(matrixFree),
    d2Slength = 0;
  else
    d2Slength = @(x,omega,m) 0;
  end  
end

% 
% the second and third summands (Sarea and Svol) are nonlinear and time consuming
% because of the partition of each voxel. In the MEX version both parts are 
% now treated in one loop.
%
% Sarea(uc) = int phiArea(A(xc+uc) / ARef)-1) dx , 
% A=surface area of triangles
if matrixFree
        if dim==3
            ARef = geometry(voxel,mOne,'A','matrixFree',matrixFree);
        else
            ARef = 0;
            alphaArea = 0;
        end
        alphas = [alphaLength alphaArea alphaVolume];
        VRef = geometry(voxel,mOne,'V','matrixFree',matrixFree);
    if alphaArea~=0 || alphaVolume~=0
        try
            [Snl, dSnl] = geometryMexC(yc,m,'S',ARef,VRef,alphas,h,doDerivative);
        catch err
            FAIRerror(err);
        end
    else
        Snl  = 0;
        dSnl = 0;
    end
else
    if alphaArea~=0 && dim==3,
        ARef = geometry(voxel,mOne,'A','matrixFree',matrixFree);
        [A,dA] = geometry(yc,m,'A','doDerivative',doDerivative,'matrixFree',matrixFree);
        ARef = kron(ARef,ones(prod(m),1));
        A = A ./ ARef;
        phiAreaFctn = str2func(phiArea);
        [H,dH,d2H] = phiAreaFctn(A,doDerivative);
        Sarea   =  (alphaArea*hd) * sum(H);
    else
        Sarea   = 0;
    end
    if alphaVolume~=0
        VRef = geometry(voxel,mOne,'V','matrixFree',matrixFree);
        [V,dV] = geometry(yc,m,'V','doDerivative',doDerivative,'matrixFree',matrixFree);
        VRef = kron(VRef,ones(prod(m),1));
        V = V ./ VRef;
        [G,dG,d2G] = phiVolume(V,doDerivative);
        Svolume   =  (alphaVolume*hd) *sum(G);
    else
        Svolume = 0;
    end
    Snl = Sarea + Svolume;
end

S = Slength + Snl;

if ~doDerivative, return; end

%------------------------------------------------------------------------------
% Compute derivatives

% Initialize derivatives
dS = reshape(dSlength,1,[]);

if matrixFree,
    d2S.solver = 'PCG-hyperElastic';
    d2S.yc     = yc;
    d2S.omega  = omega;
    d2S.m      = m;
    if alphaArea~=0 || alphaVolume~=0,
        d2Snl = @(x,m,yc) geometryMexC(yc,m,'d2S',x,ARef,VRef,alphas,h);
        d2S.diag = @(y) getDiag(omega,m,alphas,y,ARef,VRef,phiArea);
    else
        d2Snl = @(x,m,yc) 0;
        d2S.diag = @(y) 0;
    end
    dS = dS + dSnl;
    d2S.d2S = @(x,omega,m,yc) d2Slength(x,omega,m) + d2Snl(x,m,yc) ;
else
    d2S = d2Slength;
    if alphaArea~=0 && dim==3,
        dA = spdiags(1./ARef,0,length(ARef),length(ARef)) * dA;
        dSarea  = dH' * dA;
        d2Sarea = dA'* spdiags(d2H,0,length(d2H),length(d2H)) * dA;
        d2S = d2S + (hd*alphaArea)* d2Sarea;
        dS  =  dS + (hd*alphaArea) * dSarea;
    end
    if alphaVolume~=0
        dV   = spdiags(1./VRef,0,length(VRef),length(VRef)) * dV;
        %     Svolume = alpha*hd*sum(G(V(y))
        % => dSvolume = alpha*hd*e'*diag(dG)*dV = alpha*hd*dG'*dV
        dS  =  dS  + (hd*alphaVolume) * (dG' * dV);
        d2S =  d2S + (hd*alphaVolume) * (dV' * sdiag(d2G) * dV);
    end
    
end

%------------------------------------------------------------------------------
% shortcut for sparse diagonal matrices
function a = sdiag(a)
a = spdiags(reshape(a,[],1),0,length(a),length(a));

%------------------------------------------------------------------------------
% returns diagonal of d2S (interesting in matrixFree mode)
function D = getDiag(omega,m,alphas,y,ARef,VRef,phiArea)
dim   = length(omega)/2;
one   = @(i,j) One(omega,m,i,j);
h = (omega(2:2:end)-omega(1:2:end))./m;
hd    = prod((omega(2:2:end)-omega(1:2:end))./m);
if dim == 2,
  Dlin   = [reshape( one(1,1) + one(1,2),[],1)
    reshape( one(2,1) + one(2,2),[],1)];
else
  Dlin   = [ ...
    reshape(one(1,1)+(one(1,2)+one(1,3)),[],1)
    reshape(one(2,2)+(one(2,1)+one(2,3)),[],1)
    reshape(one(3,3)+(one(3,1)+one(3,2)),[],1) ];
end;
try
    D = alphas(1)*hd*Dlin + geometryMexC(y,m,'d2Sdiag',ARef,VRef,alphas,h);
catch err
    FAIRerror(err);
end

%------------------------------------------------------------------------------
function [G dG d2G] = phiVolume(x,doDerivative)
%
% phi(x) = ((x-1)^2/x)^2
%
% phi satisfies the three important conditions
%      phi(x) > 0, forall x
%      phi(1) = 0
%      phi is convex
%      phi(x) = phi(1/x)
%      phi yields det(Dy) in L_2
dG = [];
d2G = [];
% G = (x-1).*(x-1) ./x;
G = (x-1).*(x-1) ./x;
G = G.*G;
if doDerivative,
%   dG  = 1- 1./(x.*x); 
%   d2G = 2 ./ (x.*x.*x);
   dG  = 2* (x-1).^3 .* (x+1)./ x.^3; 
   d2G = 2* (x.^4-4*x+3) ./ x.^4;
end

%------------------------------------------------------------------------------
function [G dG d2G] = phiAreaConvex(x,doDerivative)
%
%   phiAreaConvex is a convex penalty function for surface area
%   regularization. This function is needed in the existence proof. In
%   order to be convex, only area growth can be penalized
%   
%
%    phi(x>=1) = 0.5 * (A/ ARef - 1)^2
%    phi(x<1 ) = 0
%
%    A    - area after deformation (24 Triangles per voxel, scalar)
%    ARef - area of reference configuration (24 Triangles per voxel, scalar)
%
dG = [];
d2G = [];
G = 0.5* (x-1).*(x-1);
G(x<1)=0;
if doDerivative,
   dG = (x-1);
   dG(x<1)=0;
   d2G = ones(size(dG));
   d2G(x<1)=0;
end

%------------------------------------------------------------------------------
function [G dG d2G] = phiAreaDoubleWell(x,doDerivative)
%
%   phiAreaDoubleWell is a penalty function for surface area
%   regularization that penalizes growth and shrinkage of area. However, 
%   this function is a double well and thus not convex. 
%   
%    phi(A) = 0.5 * (A/ ARef - 1)^2
%
%    A    - area after deformation (24 Triangles per voxel, scalar)
%    ARef - area of reference configuration (24 Triangles per voxel, scalar)
dG = [];
d2G = [];
G = 0.5* (x-1).*(x-1);
if doDerivative,
   dG = (x-1);
   d2G = ones(size(dG));
end

% matrixFree implementation of \nabla *  y and y' * \nabla
function By = lengthOperator(yc,omega,m,flag)
dim = length(omega)/2; 
h   = (omega(2:2:end)-omega(1:2:end))./m;
nn  = prod(m+1); % nodal

flag = sprintf('%s-%dD',flag,dim);
switch flag,
    case 'By-2D'
        y1  = reshape(yc(1:nn),m+1);
        y2  = reshape(yc(nn+(1:nn)),m+1);
        
        p1 = m(1) * (m(2)+1);             %    1:p1   d11 y1
        p2 = p1+(m(1)+1)*(m(2));          % p1+1:p2   d21 y1
        p3 = p2+(m(1))*(m(2)+1);          % p2+1:p3   d12 y2
        p4 = p3+(m(1)+1)*(m(2));          % p3+1:p4   d22 y2
        
        d1 = @(Y) reshape(Y(2:end,:)-Y(1:end-1,:),[],1)/h(1);
        d2 = @(Y) reshape(Y(:,2:end)-Y(:,1:end-1),[],1)/h(2);
        By = zeros(p4,1);
        By(1:p4)    = [d1(y1);d2(y1);d1(y2);d2(y2)];
        
    case 'BTy-2D'
        % the input yc is the result of B*y, therefore it can be splitted
        % into 4 parts: d11u1, d21u1,d12u2,d22u2
        
        p1 = m(1) * (m(2)+1);             %    1:p1   d11 y1
        p2 = p1+(m(1)+1)*(m(2));          % p1+1:p2   d21 y1
        p3 = p2+(m(1))*(m(2)+1);          % p2+1:p3   d12 y2
        p4 = p3+(m(1)+1)*(m(2));          % p3+1:p4   d22 y2
        
        By  = zeros(sum(2*nn),1);
        
        d1T = @(Y) reshape([-Y(1,:);Y(1:end-1,:)-Y(2:end,:);Y(end,:)],[],1)/h(1);
        d2T = @(Y) reshape([-Y(:,1),Y(:,1:end-1)-Y(:,2:end),Y(:,end)],[],1)/h(2);
        
        By(1:nn)     = d1T(reshape(yc(   1:p1),m(1),m(2)+1)) ...
            +               d2T(reshape(yc(p1+1:p2),m(1)+1,m(2)));
        
        By(nn+(1:nn)) = d1T(reshape(yc(p2+1:p3),m(1),m(2)+1)) ...
            +               d2T(reshape(yc(p3+1:p4),m(1)+1,m(2)));
        
    case 'By-3D'
        % for the linear part we need to compute the derivative of the
        % nodal grid u
        y1  = reshape(yc(1:nn),m+1);
        y2  = reshape(yc(nn+(1:nn)),m+1);
        y3  = reshape(yc(2*nn+(1:nn)),m+1);
        
        p1  = m(1)* (m(2)+1) * (m(3)+1);          %    1:p1   d11 y1
        p2  = p1+(m(1)+1)* (m(2)) * (m(3)+1);     % p1+1:p2   d21 y1
        p3  = p2+(m(1)+1)* (m(2)+1) * (m(3));     % p2+1:p3   d31 y1
        p4  = p3+(m(1))* (m(2)+1) * (m(3)+1);     % p3+1:p4   d12 y2
        p5  = p4++(m(1)+1)* (m(2)) * (m(3)+1);    % p4+1:p5   d22 y2
        p6  = p5++(m(1)+1)* (m(2)+1) * (m(3));    % p5+1:p6   d32 y2
        p7  = p6+(m(1))* (m(2)+1) * (m(3)+1);     % p6+1:p7   d13 y3
        p8  = p7+(m(1)+1)* (m(2)) * (m(3)+1);     % p7+1:p8   d23 y3
        p9  = p8+(m(1)+1)* (m(2)+1) * (m(3));     % p8+1:p9   d33 y3
        By = zeros(p9,1);
        
        d1 = @(Y) reshape(Y(2:end,:,:)-Y(1:end-1,:,:),[],1)/h(1);
        d2 = @(Y) reshape(Y(:,2:end,:)-Y(:,1:end-1,:),[],1)/h(2);
        d3 = @(Y) reshape(Y(:,:,2:end)-Y(:,:,1:end-1),[],1)/h(3);
        
        By(1:p9)    = [...
            d1(y1);d2(y1);d3(y1);...
            d1(y2);d2(y2);d3(y2);...
            d1(y3);d2(y3);d3(y3) ];
        
    case 'BTy-3D'
        % the input yc is the result of B*y, therefore it can be splitted into 10
        % parts: d11y1, d21y1,d31y1,d12y2,d22y2,d32y2,d13y3,d23y3,u33y3
        
        p1  = m(1)* (m(2)+1) * (m(3)+1);          %    1:p1   d11 y1
        p2  = p1+(m(1)+1)* (m(2)) * (m(3)+1);     % p1+1:p2   d21 y1
        p3  = p2+(m(1)+1)* (m(2)+1) * (m(3));     % p2+1:p3   d31 y1
        p4  = p3+(m(1))* (m(2)+1) * (m(3)+1);     % p3+1:p4   d12 y2
        p5  = p4+(m(1)+1)* (m(2)) * (m(3)+1);     % p4+1:p5   d22 y2
        p6  = p5+(m(1)+1)* (m(2)+1) * (m(3));     % p5+1:p6   d32 y2
        p7  = p6+(m(1))* (m(2)+1) * (m(3)+1);     % p6+1:p7   d13 y3
        p8  = p7+(m(1)+1)* (m(2)) * (m(3)+1);     % p7+1:p8   d23 y3
        p9  = p8+(m(1)+1)* (m(2)+1) * (m(3));     % p8+1:p9   d33 y3
        By  = zeros(3*nn,1);
        
        d1T = @(Y) reshape(d1t(Y),[],1)/h(1);
        d2T = @(Y) reshape(d2t(Y),[],1)/h(2);
        d3T = @(Y) reshape(d3t(Y),[],1)/h(3);
        
        
        By(1:nn)   = d1T(reshape(yc(   1:p1),m(1),(m(2)+1),(m(3)+1))) ...
            +             d2T(reshape(yc(p1+1:p2),m(1)+1,m(2),m(3)+1)) ...
            +             d3T(reshape(yc(p2+1:p3),(m(1)+1),(m(2)+1),(m(3))));
        
        By(nn+(1:nn)) = ...
            d1T(reshape(yc(p3+1:p4),m(1),m(2)+1,m(3)+1)) ...
            +   d2T(reshape(yc(p4+1:p5),m(1)+1,m(2),m(3)+1)) ...
            +   d3T(reshape(yc(p5+1:p6),m(1)+1,m(2)+1,m(3)));
        
        By(2*nn+(1:nn)) = ...
            d1T(reshape(yc(p6+1:p7),m(1),m(2)+1,m(3)+1)) ...
            +   d2T(reshape(yc(p7+1:p8),m(1)+1,m(2),m(3)+1)) ...
            +   d3T(reshape(yc(p8+1:p9),m(1)+1,m(2)+1,m(3)));
    otherwise
        error('Unknown flag %s', flag);
end

% partial derivatives with respect to x_1 (mf)
function y = d1t(Y)
m = size(Y);
y = zeros(m+[1,0,0]);
y(1,:,:) = -Y(1,:,:);
y(2:end-1,:,:) = Y(1:end-1,:,:)-Y(2:end,:,:);
y(end,:,:) = Y(end,:,:);

% partial derivatives with respect to x_2 (mf)
function y = d2t(Y)
m = size(Y);
y = zeros(m+[0,1,0]);
y(:,1,:) = -Y(:,1,:);
y(:,2:end-1,:) = Y(:,1:end-1,:)-Y(:,2:end,:);
y(:,end,:) = Y(:,end,:);

% partial derivatives with respect to x_3 (mf)
function y = d3t(Y) 
m = size(Y); if length(m) == 2, m = [m,1]; end;
y = zeros(m+[0,0,1]);
y(:,:,1) = -Y(:,:,1);
y(:,:,2:end-1) = Y(:,:,1:end-1)-Y(:,:,2:end);
y(:,:,end) = Y(:,:,end);

% for diagonal of legnth operator
function o = One(omega,m,i,j)
h = (omega(2:2:end)-omega(1:2:end))./m;
m = m + 1;
o = ones(m)/h(j)^2;
switch j,
  case 1, o(2:end-1,:,:) = 2*o(2:end-1,:,:);
  case 2, o(:,2:end-1,:) = 2*o(:,2:end-1,:);
  case 3, o(:,:,2:end-1) = 2*o(:,:,2:end-1);
end;

%------------------------------------------------------------------------------
function runMinimalExample
help(mfilename);
% specify elasticity parameters
alphaLength = 1.4;
alphaArea   = 4.3;
alphaVolume = 5.7;

% --- 2D --- 
omega = [0 2 0 2];
m     = [17 19];
xc    = getNodalGrid(omega,m);
uc    = 1e-2*randn(size(xc));

Hyper  = @(uc,omega,m,matrixFree) feval(mfilename,uc,omega,m,...
'alphaLength',alphaLength,'alphaArea',alphaArea,'alphaVolume',alphaVolume,...
'matrixFree',matrixFree);


[mbS, mbdS, mbd2S] = Hyper(uc,omega,m,0);
[mfS, mfdS, mfd2S] = Hyper(uc,omega,m,1);

fprintf('HyperElastic potential : %f (mb) | %f (mf) \n' , mbS,mfS);
RE = norm(mbdS-mfdS)/norm(mbdS);
fprintf('rel.error dS (mb vs. mf) : %e \n' , RE);
x = randn(2*prod(m+1),1);
RE = norm(mbd2S*x-mfd2S.d2S(x,omega,m,xc+uc))/norm(mbd2S*x);
fprintf('rel.error d2S (mb vs. mf) : %e \n' , RE);
RE = norm(full(diag(mbd2S))-mfd2S.diag(xc+uc))/norm(full(diag(mbd2S)));
fprintf('rel.error diag (mb vs. mf) : %e \n' , RE);


figure(1)
subplot(1,2,1);
spy(mbd2S)
title(sprintf('%s-spy(d2S)-(2D)',mfilename));  


% --- 3D --- 
omega = [0 2 0 3 0 2];
m     = [6 5 7];
xc    = getNodalGrid(omega,m);
uc    = 1e-2*randn(size(xc));
[mbS, mbdS, mbd2S] = Hyper(uc,omega,m,0);
[mfS, mfdS, mfd2S] = Hyper(uc,omega,m,1);

  
fprintf('HyperElastic potential : %f (mb) | %f (mf) \n' , mbS,mfS);
RE = norm(mbdS-mfdS)/norm(mbdS);
fprintf('rel. error dS (mb vs. mf) : %e \n' , RE);
x = randn(3*prod(m+1),1);
RE = norm(mbd2S*x-mfd2S.d2S(x,omega,m,xc+uc))/norm(mbd2S*x);
fprintf('rel.error d2S (mb vs. mf) : %e \n' , RE);
RE = norm(full(diag(mbd2S))-mfd2S.diag(xc+uc))/norm(full(diag(mbd2S)));
fprintf('rel.error diag (mb vs. mf) : %e \n' , RE);

figure(1)
subplot(1,2,2);
spy(mbd2S)
title(sprintf('%s-spy(d2S)-(3D)',mfilename));
%==============================================================================
  

