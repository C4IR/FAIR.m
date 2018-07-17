%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function [Y,LM] = LMreg(type,LM,X,varargin)
%
% Landmark based Image Registration
%
% Input:
%   type         in {'linear','quadratic','TPS'} denotes the transformation model
%   LM          the location of landmarks m-by-2*d, d = dimension
%               t = LM(:,1:d), LM in template, r = LM(:,d+1:2*d) LM in reference
%   X           grid where the LM solution is evaluated
%   varargin    optional arguments, like theta for TPS, see below
%
% Output:
%    Y                 transformed grid
%    LM         updated landmarks
%
% see also E5_Hands_TPS for an example
%==============================================================================

function [Y,LM] = LMreg(type,LM,X,varargin)

if nargin==0
    help(mfilename)
    runMinimalExample; 
    return;
end

theta = 0;
for k=1:2:length(varargin), % overwrites default parameter 
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

% reshape grid and extract landmarks for T and R
dim = size(LM,2)/2;
t   = LM(:,1:dim);
r   = LM(:,dim+1:2*dim);
X   = reshape(X,[],dim);
Y   = zeros(size(X));

switch type,
  case 'linear',
    % linear transformation model
    Q  = [r,ones(size(LM,1),1)];
    w  = (Q'*Q)\(Q'*t);
    for j = 1:dim,
      Y(:,j)        = [X,ones(size(X,1),1)]*w(:,j);
      LM(:,2*dim+j) = Q*w(:,j);      
    end;
  case 'quadratic',
    % quadratic 2D transformation model
    Q = [ones(size(LM,1),1),LM(:,[3:4]),LM(:,3).^2,LM(:,4).^2,LM(:,3).*LM(:,4)];
    w = (Q'*Q)\(Q'*LM(:,1:2));
    m = @(w,X) (w(1)+w(2)*X(:,1)+w(3)*X(:,2)...
      +w(4)*X(:,1).^2+w(5)*X(:,2).^2+w(6)*X(:,1).*X(:,2)); 
    for j = 1:dim,
      Y(:,j)        = [ones(size(X,1),1),X,X(:,1).^2,X(:,2).^2,X(:,1).*X(:,2)]*w(:,j);
      LM(:,2*dim+j) = Q*w(:,j);      
    end;
  case 'TPS',
    % Thin-Plate-Spline transformation model, 
    w = getTPScoefficients(LM(:,1:4),'theta',theta);
    [Y,yLM] = evalTPS(LM,w,X);
    LM(:,[5,6]) = reshape(yLM,[],2);
  otherwise
    fprintf('unknown model %s in <%s>\n',type,mfilename);
end;


if (nargout < 2) || (dim ~= 2), return; end;

% approximate the inverse transformation on the landmarks
Y = reshape(Y,[],2);
LM(:,7) = griddata(Y(:,1),Y(:,2),X(:,1),LM(:,1),LM(:,2));
LM(:,8) = griddata(Y(:,1),Y(:,2),X(:,2),LM(:,1),LM(:,2));
Y = Y(:);

%------------------------------------------------------------------------------

function runMinimalExample
setup2DhandData
omegaT = omega(1,:);
omegaR = omega(end,:);
xT = getCellCenteredGrid(omegaT,m);
xR = getCellCenteredGrid(omegaR,m);
Tc = imgModel(dataT,omegaT,xT);
Rc = imgModel(dataR,omegaR,xR);

theta = 10;
[yc,LM] = LMreg('TPS',LM(:,1:4),xR,'theta',theta);
TLM = imgModel(dataT,omegaT,yc);


FAIRfigure(1,'figname',mfilename); clf;

subplot(1,3,1); 
viewImage(Tc,omegaT,m); hold on; axis off;
ph = plotLM(LM(:,1:2),'numbering','on','color','r');
set(ph,'linewidth',2,'markersize',20);
title(sprintf('%s','T&LM'),'fontsize',20);

subplot(1,3,2); 
viewImage(Rc,omegaR,m); hold on;  axis off;
ph = plotLM(LM(:,3:4),'numbering','on','color','g','marker','+');
set(ph,'linewidth',2,'markersize',20);
title(sprintf('%s','R&LM'),'fontsize',20);

subplot(1,3,3); 
viewImage(TLM,omegaR,m); hold on;  axis off;
ph = plotLM(LM(:,3:4),'numbering','off','color','g','marker','+');
qh = plotLM(LM(:,7:8),'numbering','off','color','m','marker','x');
rh = plot(LM(:,[3,7])',LM(:,[4,8])','m-','linewidth',3);
set([ph;qh;rh],'linewidth',2,'markersize',20);
title(sprintf('T(Y^{TPS},\\theta=%s)&LM',num2str(theta)),'fontsize',30);

%==============================================================================
