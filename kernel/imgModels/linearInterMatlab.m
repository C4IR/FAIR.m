%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% 
% function [Tc,dT] = linearInterMatlab(T,omega,x,varargin)
%
% MATLAB's linear interpolator 
% for the data T given on a cell-centered grid evaluated at x
% uses finite difference approximation for derivatives
% see Exercise 3.4
%==============================================================================

function [Tc,dT] = linearInterMatlab(T,omega,x,varargin)
         
Tc = mfilename('fullpath'); dT = []; 

if nargin == 0, 
  help(mfilename);
  testOneImgModel(mfilename);
  return;
elseif nargin == 1 && isempty(T),
  return;
end;
      
% if nargin == 0, return filename
if nargin == 0, 
  Tc = mfilename('fullpath'); 
  dT = []; 
  return; 
end;


% flag for computing the derivative
doDerivative = (nargout>1);

for k=1:2:length(varargin), % overwrite default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

% get data size m, cell size h, dimension dim, and number n of interpolation points
dim = length(omega)/2;
m   = size(T);         if dim == 1, m = numel(T); end;
h   = (omega(2:2:end)-omega(1:2:end))./m; 
n   = numel(x)/dim;
x   = reshape(x,n,dim);
Tc  = zeros(n,1);
dT  = sparse(n,dim*n);

xi = @(i) (omega(2*i-1)+h(i)/2:h(i):omega(2*i)-h(i)/2)';
switch dim
    case 1
        xc = getCellCenteredGrid(omega,m);
        Tc = reshape(interp1(xc,T,x,'linear'),[],1);         Tc(isnan(Tc)) = 0;
        if ~doDerivative, return; end;
        h = (omega(2:2:end)-omega(1:2:end))./(length(T)*10);
        Tp = interp1(xc,T,x+h(1)/2,'linear');  Tp(isnan(Tp)) = 0;
        Tm = interp1(xc,T,x-h(1)/2,'linear');  Tm(isnan(Tm)) = 0;
        dT = (Tp-Tm)/h(1);
    case 2
        [X1,X2] = meshgrid(xi(1),xi(2));
        Tfctn = @(x1,x2) reshape(interp2(X1,X2,permute(T,[2,1]),x1,x2,'linear'),[],1);
        Tc = Tfctn(x(1:n),x(n+1:end));        Tc(isnan(Tc)) = 0; 
        if ~doDerivative, return; end;
        Tp = Tfctn(x(1:n)+h(1)/2,x(n+1:end)); Tp(isnan(Tp)) = 0;
        Tm = Tfctn(x(1:n)-h(1)/2,x(n+1:end)); Tm(isnan(Tm)) = 0;
        dT(:,1) = (Tp-Tm)/h(1);
        Tp = Tfctn(x(1:n),x(n+1:end)+h(2)/2); Tp(isnan(Tp)) = 0;
        Tm = Tfctn(x(1:n),x(n+1:end)-h(2)/2); Tm(isnan(Tm)) = 0;
        dT(:,2) = (Tp-Tm)/h(2);
    case 3
        [X1,X2,X3] = meshgrid(xi(1),xi(2),xi(3));
        Tfctn = @(x1,x2,x3) reshape(...
          interp3(X1,X2,X3,double(permute(T,[2,1,3])),x1,x2,x3,'linear'),[],1);
        Tc = Tfctn(x(1:n),x(n+1:2*n),x(2*n+1:end));         Tc(isnan(Tc)) = 0;
        if ~doDerivative, return; end;
        Tp = Tfctn(x(1:n)+h(1)/2,x(n+1:2*n),x(2*n+1:end));  Tp(isnan(Tp)) = 0;
        Tm = Tfctn(x(1:n)-h(1)/2,x(n+1:2*n),x(2*n+1:end));  Tm(isnan(Tm)) = 0;
        dT(:,1) = (Tp-Tm)/h(1);
        Tp = Tfctn(x(1:n),x(n+1:2*n)+h(2)/2,x(2*n+1:end));  Tp(isnan(Tp)) = 0;
        Tm = Tfctn(x(1:n),x(n+1:2*n)-h(2)/2,x(2*n+1:end));  Tm(isnan(Tm)) = 0;
        dT(:,2) = (Tp-Tm)/h(2);
        Tp = Tfctn(x(1:n),x(n+1:2*n),x(2*n+1:end)+h(3)/2);  Tp(isnan(Tp)) = 0;
        Tm = Tfctn(x(1:n),x(n+1:2*n),x(2*n+1:end)-h(3)/2);  Tm(isnan(Tm)) = 0;
        dT(:,3) = (Tp-Tm)/h(3);
end
if doDerivative,
    dT = spdiags(dT,n*(0:(dim-1)),n,dim*n);
end
%==============================================================================
