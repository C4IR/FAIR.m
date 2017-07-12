%==============================================================================
% This code is part of the Matlab-based toolbox 
% LagLDDDM - A Lagrangian Gauss--Newton--Krylov Solver for Mass- and 
%                        Intensity-Preserving Diffeomorphic Image Registration
% 
% For details and license info see 
% - https://github.com/C4IR/FAIR.m/tree/master/add-ons/LagLDDMM
%
% function [yx,Dyx] = linearInterGrid(yc,omega,m,x)
%
% interpolates grid yc at positions x
%
% Input:
% 
% yc    - grid
% omega - spatial domain
% m     - number of discretization points
% x     - interpolation points
% 
% Output: 
% 
% yx    - y(x), yc evaluated at x
% Dyx   - derivative of y(x) w.r.t. x
%==============================================================================
function [yx, Dyx] = linearInterGrid(yc,omega,m,x,varargin)

if nargin==0
    runMinimalExample;
    return;
end;
inter = @linearInterMex;
doDerivative = (nargout==2);
for k=1:2:length(varargin)    % overwrites default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

dim = length(omega)/2;
h   = (omega(2:2:end)-omega(1:2:end))./m;

switch numel(yc)
    case prod(m)*dim % cell-centered
        yc = reshape(yc,[],dim);
        yx = reshape(zeros(size(x)),[],dim);
        Dyx = [];
        for d=1:dim
            [yx(:,d),Dyt] = inter(reshape(yc(:,d),m),omega,x(:),'doDerivative',doDerivative);
            if doDerivative
                Dyx = [Dyx ; Dyt];
            end
        end
    case prod(m+1)*dim % nodal
        yc = reshape(yc,[],dim);
        omega = omega + reshape([-h;h]/2,1,[]);
        yx = reshape(zeros(size(x)),[],dim);
        Dyx = [];
        for d=1:dim
            if dim==1
                ycoeff = yc(:);
            else
                ycoeff = reshape(yc(:,d),m+1);
            end
            [yx(:,d),Dyt] = inter(ycoeff,omega,x(:),'doDerivative',doDerivative);
            if doDerivative
                Dyx = [Dyx ; Dyt];
            end
        end
    case sum(prod(ones(dim,1)*m+eye(dim),2)) % staggered
        yx = reshape(zeros(size(x)),[],dim);
        ns  = prod(ones(length(m),1)*m+eye(length(m)),2); %staggered
        pd  = 0;
        Dyx = [];
        for d=1:dim
            omegaGrid = omega;
            omegaGrid(2*d-1:2*d) = omegaGrid(2*d-1:2*d) + [-h(d), h(d)]/2;
            mGrid = m + ((1:dim)==d);
            [yx(:,d),Dyt] = inter(reshape(yc(pd+ (1:ns(d))),mGrid),omegaGrid,x(:),'doDerivative',doDerivative);
            if doDerivative
                Dyx = [Dyx ; Dyt];
            end
            pd = pd+ns(d);
        end
    otherwise
        error('cannot deal this grid!');
end
yx = yx(:);
          


function runMinimalExample

omega = [0,10,0,8]; m = [8,9]; p = [5,6];
w  = zeros([p,2]);  w(3,3,1) = 0.07; w(3,4,2) = -0.06;
xc = getNodalGrid(omega,m);
yc  = splineTransformation2D(w(:),xc,'omega',omega,'m',m+1,'p',p,'Q',[]);

% some points describing a rectangle
X = [5;6;5;6;2;2;4;4;];

yX = feval(mfilename,yc,omega,m,X);

fctn = @(x) linearInterGrid(yc,omega,m,x);
checkDerivative(fctn,X(:)+randn(size(X(:)))/10);

figure(1); clf;
subplot(1,2,1);
plotGrid(xc,omega,m); hold on; 
plot(X(1:4),X(5:end),'rx');
title('initial grid + points');
subplot(1,2,2);
plotGrid(yc,omega,m); hold on; 
plot(yX(1:4),yX(5:end),'rx');
title('transformed grid + transformed points');