%==============================================================================
% ##1
%
% function [T,dT] = getPICMatrixAnalyticIntegral(omega,m,hp,xp,eps,varargin)
%
% Building push-forward matrix for aParticle-In-Cell method using anlytic
% integration.
%
%
% Input:
%  omega    - representation of computational domain
%  mc       - discretization size of sampling grid of data
%  mf       - discretization size of high-resolution image
%  hp       - cell size on particle grid used for computing particle's mass
%             which is rho_i \approx prod(hp)* \rho(x_i)
%  xp       - particle positions, size(xp) = [dim*np,1];
%  eps      - width of particles
%  rho      - particles mass, size(rho)=[np,1]
%  varargin - optional additional input
%
% Output:
%  T        - push-forward matrix, i.e. rho(xp) = C*rho
%  dT       - derivative of (C(xp)*rho) with respect to xp.
%
% =========================================================================
function [T,dT] = getPICMatrixAnalyticIntegral(omega,mc,mf,xp,varargin)

if nargin==0, help(mfilename); runMinimalExample; return; end

hp = (omega(2:2:end)-omega(1:2:end))./mf;   % cell-size in particle mesh
epsP = hp;                                  % width of particles
doDerivative = (nargout==2);
for k=1:2:length(varargin) % overwrites defaults
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

dT = [];
dim = numel(omega)/2;
xp  = reshape(xp,[],dim);
np  = size(xp,1);         % number of particles
n   = prod(mc);            % number of voxels in sampling grid

h   = (omega(2:2:end)-omega(1:2:end))./mc;
pwidth = (ceil(epsP./h));    % upper bound for support of basis functions

% map particle positions to the domain [0,m(1)] x [0,m(dim)]
for i=1:dim
    xp(:,i) = (xp(:,i)-omega(2*i-1))/h(i);
end
% omega = zeros(1,2*dim);omega(2:2:end) = mc;

% get cell index of particles center of mass
P = ceil(xp);
w = xp-(P-1);

switch dim
    case 1
        B = reshape(int1D(w,pwidth,epsP,h,hp),[],1);
        J = repmat(  (1:np)',[2*pwidth+1,1]);
        I = reshape(bsxfun(@plus,P(:,1),-pwidth:pwidth),[],1);
        
        
        valid = (1<=I) & (I<=mc(1));
        I = I(valid); J = J(valid); B = B(valid);
        
        T = sparse(I,J,B,n,np);
        if doDerivative
            bx = reshape(diff1D(w,pwidth,epsP,h,hp),[],1);
            ids = (1:numel(bx)); ids = ids(valid);
            dT = @(rho) sparse(I,J,rho(J).*bx(ids),n,np)/h;
        end
        
    case 2
        s2i  = @(i) i(:,1)+ mc(1)   *(i(:,2)-1);   % sub2ind for cell-centered grid
        
        B1 = int1D(w(:,1),pwidth(1),epsP(1),h(1),hp(1));
        B2 = int1D(w(:,2),pwidth(2),epsP(2),h(2),hp(2));
        if doDerivative
            b1 = diff1D(w(:,1),pwidth(1),epsP(1),h(1),hp(1));
            b2 = diff1D(w(:,2),pwidth(2),epsP(2),h(2),hp(2));
        end
        
        
        nVoxel = prod([size(B1,2),size(B2,2)]);
        J = repmat((1:np)',[nVoxel 1]);
        I = zeros([np*nVoxel,2]); B = zeros(np*nVoxel,1); pp = 1;
        bx = B; by = B;
        for px = -pwidth(1):pwidth(1)
            for py = -pwidth(2):pwidth(2)
                idx = (pp-1)*np+(1:np); pp = pp+1;
                I(idx,:) = [P(:,1)+px, P(:,2)+py];
                
                B(idx) =  B1(:,px+pwidth(1)+1).*B2(:,py+pwidth(2)+1);
                
                if doDerivative
                    bx(idx) = B2(:,py+pwidth(2)+1).*b1(:,px+pwidth(1)+1);
                    by(idx) = B1(:,px+pwidth(1)+1).*b2(:,py+pwidth(2)+1);
               end
            end
        end
        valid = (1<=I(:,1)) & (I(:,1)<=mc(1)) & (1<=I(:,2)) & (I(:,2)<=mc(2));
        I = s2i(I(valid,:));
        J = J(valid);
        B = B(valid);
        T = sparse(I,J,B,n,np);
        
        if doDerivative
            ids = (1:numel(bx))'; ids = ids(valid);
            dT = @(rho) [ sparse(I,J,rho(J).*bx(ids)/h(1),n,np), sparse(I,J,rho(J).*by(ids)/h(2),n,np)];
        end
    case 3
        % sub2ind for cell-centered grid
        s2i  = @(i) i(:,1)+ mc(1)   *(i(:,2)-1)+(mc(1)  )*(mc(2)  )*(i(:,3)-1);
        
        B1 = int1D(w(:,1),pwidth(1),epsP(1),h(1),hp(1));
        B2 = int1D(w(:,2),pwidth(2),epsP(2),h(2),hp(2));
        B3 = int1D(w(:,3),pwidth(3),epsP(3),h(3),hp(3));
        
        
        nVoxel = prod([size(B1,2),size(B2,2),size(B3,2)]);
        I = zeros([np*nVoxel,3]); J = repmat( (1:np)',[nVoxel,1]); B = zeros(np*nVoxel,1);
        pp = 1;
        if doDerivative
            b1 = diff1D(w(:,1),pwidth(1),epsP(1),h(1),hp(1));
            b2 = diff1D(w(:,2),pwidth(2),epsP(2),h(2),hp(2));
            b3 = diff1D(w(:,3),pwidth(3),epsP(3),h(3),hp(3));
            
            bx = B; by = B; bz = B;
        end
        for pz = -pwidth(3):pwidth(3)
            i3 = P(:,3)+pz;
            B3t = B3(:,pz+pwidth(3)+1);
            for py = -pwidth(2):pwidth(2)
                i2 = P(:,2)+py;
                B2t = B2(:,py+pwidth(2)+1);
                B23 = B3t.*B2t;
                for px = -pwidth(1):pwidth(1)
                    idx = (pp-1)*np+(1:np); pp = pp+1;
                    I(idx,1) = P(:,1)+px;
                    I(idx,2)=i2;
                    I(idx,3)=i3;
                    
                    % remove cells that lie out of the domain
                    B1t = B1(:,px+pwidth(1)+1);
                    
                    B(idx)   =  B1t.*B23;
                    
                    if doDerivative
                       % compute derivatives of weights
                        bx(idx)     = B23.*b1(:,px+pwidth(1)+1);
                        by(idx)     = B1t.*B3t.*b2(:,py+pwidth(2)+1);
                        bz(idx)     = B1t.*B2t.*b3(:,pz+pwidth(3)+1);
                    end
                end
            end
        end
        valid =    (1<=I(:,1)) & (I(:,1)<=mc(1)) ...
            & (1<=I(:,2)) & (I(:,2)<=mc(2)) ...
            & (1<=I(:,3)) & (I(:,3)<=mc(3));
        
        I = s2i(I(valid,:));
        J = J(valid);
        B = B(valid);
        T = sparse(I,J,B,n,np);
        
        if doDerivative
            ids = (1:numel(bx))'; ids = ids(valid);
            D1 = @(rho) sparse(I,J,rho(J).*bx(ids)/h(1),n,np);
            D2 = @(rho) sparse(I,J,rho(J).*by(ids)/h(2),n,np);
            D3 = @(rho) sparse(I,J,rho(J).*bz(ids)/h(3),n,np);
            dT = @(rho) [D1(rho), D2(rho), D3(rho)];
        end
        
    otherwise
        error('dimension must be 1,2,3')
end

function Bij = int1D(w,pwidth,eps,h,hp)
Bij = zeros(numel(w),2*pwidth+1);
Bleft = B(-pwidth-w,eps,h);
for p = -pwidth:pwidth
    Bright = B(1+p-w,eps,h);
    Bij(:,p+pwidth+1)  = hp*(Bright - Bleft);
    Bleft = Bright;
end
function Bij = diff1D(w,pwidth,eps,h,hp)
Bij = zeros(numel(w),2*pwidth+1);
Bleft = b(-pwidth-w,eps,h);
for p = -pwidth:pwidth
    Bright = b(1+p-w,eps,h);
    Bij(:,p+pwidth+1)  = hp*(Bright - Bleft);
    Bleft = Bright;
end
Bij = - Bij;


function bij = b(x,eps,h)
bij = zeros(numel(x),1);

ind1 = (-eps/h<=x)&(x<=0);
ind2 = (0<x)&(x<=eps/h);

bij(ind1)  = 1 + h*x(ind1)./eps;
bij(ind2)  = 1 - h*x(ind2)./eps;
bij = bij /eps;

function Bij = B(x,eps,h)
Bij = zeros(numel(x),1);

ind1 = (-eps/h<=x)&(x<=0);
ind2 = (0<x)&(x<=eps/h);
ind3 = (eps/h<x);

Bij(ind1) = x(ind1) + 1./(2*eps/h).*x(ind1).^2+eps/(h*2);
Bij(ind2) = x(ind2) - 1./(2*eps/h).*x(ind2).^2+eps/(h*2);
Bij(ind3) = eps/h;
Bij = Bij/eps;



function [T,dT] = derivativeTestFctn(omega,m,hc,xp,rho)
[T,dT] = feval(mfilename,omega,m,hc,xp);
T = T*rho;
dT = dT(rho);

function runMinimalExample

%  ========== 1 D ==================
omega = [0 1]; mc = 64; mf = 64;
hc     = (omega(2:2:end)-omega(1:2:end))./mc;
hf     = (omega(2:2:end)-omega(1:2:end))./mf;
rho   = ones(mc,1); rho([1:8,end-8:end]) = 0;

xp = getCellCenteredGrid(omega,mc);
xp = xp +.02;

% transport rho
C = feval(mfilename,omega,mf,mc,xp);
rhonew = C*rho;

% visualize result
figure(1);clf;
subplot(1,2,1)
plot(getCellCenteredGrid(omega,mc),rho);
title(sprintf('rho, mass:%e',prod(hc)*sum(rho)));

subplot(1,2,2)
plot(getCellCenteredGrid(omega,mf),rhonew)
title(sprintf('rhonew, mass:%e',prod(hf)*sum(rhonew)));
% check derivative
fctn = @(xp) derivativeTestFctn(omega,mf,mc,xp,rho);
checkDerivative(fctn,xp(:));
%  ========== 2 D ==================
omega = [-2 2 -2 2]; mc = [32 32]; mf = mc*2;
hc     = (omega(2:2:end)-omega(1:2:end))./mc;
hf     = (omega(2:2:end)-omega(1:2:end))./mf;
f      = @(x) 0 + (sqrt(x(:,1).^2 + .4*x(:,2).^2)<.7);
rho    = f(reshape(getCellCenteredGrid(omega,mc),[],2));

% choose parameters
xp = getCellCenteredGrid(omega,mc);
xp = xp+.1*hc(1) ;

% transport rho
C = feval(mfilename,omega,mf,mc,xp);
rhonew = C*rho;

% visualize result
figure(2); clf;
subplot(1,2,1);
viewImage2Dsc(rho,omega,mc);
title(sprintf('rho, mass:%e',prod(hc)*sum(rho)));

subplot(1,2,2);
viewImage2Dsc(rhonew,omega,mf);
title(sprintf('rho, mass:%e',prod(hf)*sum(rhonew)));

% check derivative
fctn = @(xp) derivativeTestFctn(omega,mf,mc,xp,rho);
checkDerivative(fctn,xp(:));


%  ========== 3 D ==================
omega = [-2 2 -2 2 -3 3]; mc = [16 32 24]; mf = mc;
hc     = (omega(2:2:end)-omega(1:2:end))./mc;
hf     = (omega(2:2:end)-omega(1:2:end))./mf;
f      = @(x) 0 + (sqrt(x(:,1).^2 + .3*x(:,2).^2 + .2*x(:,3).^2)<.7);
rho    = f(reshape(getCellCenteredGrid(omega,mc),[],3));

% choose parameters
xp = getCellCenteredGrid(omega,mc);
xp = xp+0.1*hc(1) ;

% transport rho
C = feval(mfilename,omega,mf,mc,xp);
rhonew = C*rho;


% visualize result
figure(2); clf;
subplot(1,2,1);
imgmontage(rho,omega,mc);
title(sprintf('rho, mass:%e',prod(hc)*sum(rho)));

subplot(1,2,2);
imgmontage(rhonew,omega,mf);
title(sprintf('rho, mass:%e',prod(hf)*sum(rhonew)));

% check derivative
fctn = @(xp) derivativeTestFctn(omega,mf,mc,xp,rho);
checkDerivative(fctn,xp(:));


