%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function [Y,Z,c] = TPSmex(LM,X,varargin)
%
% Landmark Based Thin-Plate-Spline Transformation
% minimizes sum_{i} { ||Y^i(r_j)-t^i_j|| + theta*S(Y^i) } != min
% where t and r denote the landmarks in template and reference
%
% the solution is explicitly known:
%   y = sum_{j} lambda_j rho(|x-r_j|) + p1(x),
% where rho is the dimension dependent radial basis function and the coefficients 
% a of p1 and the lambda satisfy (A(j,k)=rho(|r_j,r_k|), C = [1,r])
% | A  C | * |lambda| = |t|
% | C' 0 |   |a     |   |0|
%
% Input:
%   LM       position of landmarks [t(:,1),...,t(:,d),r(:,1),...,r(:,d)]
%   X        location of points the TPS has to be evaluated
%   [theta]  optional smoothing parameter, default: theta=0
% Output
%   Y = TPS(X)
%   Z = TPS(r)
%   c = TPS coefficients
%
% cf. YTPS for examples
%==============================================================================

function [Y,Z,c] = TPSmex(LM,X,varargin)


if nargin==0
    runMinimalExample; return;
end

theta = [];
omega = [];
m     = [];
for k=1:2:length(varargin), % overwrite default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

[L,dim] = size(LM);
dim = dim/2;

if (numel(X) == 1)
    [Y,Z,c] = TPSmexC(LM(:,1:dim),LM(:,dim+1:2*dim),X,theta,omega,m);
else
    [Y,Z,c] = TPSmexC(LM(:,1:dim),LM(:,dim+1:2*dim),X,theta);
end
%------------------------------------------------------------------------------
function runMinimalExample
setup2DhandData
omegaT = omega(1,:); %#ok<NODEF>
omegaR = omega(end,:);
xT = getCellCenteredGrid(omegaT,m);
xR = getCellCenteredGrid(omegaR,m);
Tc = imgModel(dataT,omegaT,xT);
Rc = imgModel(dataR,omegaR,xR);

theta = 10;
[yc,z] = TPSmex(LM(:,1:4),xR,'theta',theta);
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
  z = reshape(z,[],2);
  viewImage(TLM,omegaR,m); hold on;  axis off;
  ph = plotLM(LM(:,3:4),'numbering','off','color','g','marker','+');
  qh = plotLM(z,'numbering','off','color','m','marker','x');
  rh = plot([LM(:,3) z(:,1)]',[LM(:,4) z(:,2)]','m-','linewidth',3);
  set([ph;qh;rh],'linewidth',2,'markersize',20);
  title(sprintf('T(Y^{TPS},\\theta=%s)&LM',num2str(theta)),'fontsize',30);

help(mfilename);
%==============================================================================
