%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function [Tc,dT] = splineInter(T,omega,x,varargin);
%
% cubic B-spline interpolator for the data T given on a cell-centered grid 
% evaluated at x, see FAIR 3.4 p26.
%
% note: T are the spline coefficients, img(x) = \sum T(j) b_j(x)
% use get getSplineCoefficients.m to compute the coefficients
%==============================================================================

function [Tc,dT] = splineInter(T,omega,x,varargin)
         
Tc = mfilename('fullpath'); dT = []; 

if nargin == 0, 
  help(mfilename)
  testOneImgModel(mfilename);
  return;
elseif nargin == 1 && isempty(T),
  return;
end;

% flag for computing the derivative
doDerivative = (nargout>1);
matrixFree   = 0;
for k=1:2:length(varargin), % overwrite default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

% get data size m, cell size h, dimension d, and number n of interpolation points
dim = length(omega)/2;
m   = size(T);         if dim == 1, m = numel(T); end;
h   = (omega(2:2:end)-omega(1:2:end))./m; 
n   = length(x)/dim;    
x   = reshape(x,n,dim);
% map x from [h/2,omega-h/2] -> [1,m],
for i=1:dim, x(:,i) = (x(:,i)-omega(2*i-1))/h(i) + 0.5; end;

Tc = zeros(n,1); dT = [];                   % initialize output
if doDerivative, dT = zeros(n,dim);  end;   % allocate memory in column format
% BIG-fix 2011
Valid = @(j) (-1<x(:,j) & x(:,j)<m(j)+2);   % determine indices of valid points

switch dim,
  case 1, valid = find( Valid(1) );   
  case 2, valid = find( Valid(1) & Valid(2) );   
  case 3, valid = find( Valid(1) & Valid(2) & Valid(3) );   
end;

if isempty(valid),                        
  if doDerivative, 
      dT = sparse(n,dim*n); 
  end; % allocate memory incolumn format
  return; 
end;

pad = 3; TP = zeros(m+2*pad);             % pad data to reduce cases

P = floor(x); x = x-P;                    % split x into integer/remainder
p = @(j) P(valid,j); xi = @(j) x(valid,j);
% increments for linearized ordering
i1 = 1; i2 = size(T,1)+2*pad; i3 = (size(T,1)+2*pad)*(size(T,2)+2*pad);

b0  = @(j,xi) motherSpline(j,xi);         % shortcuts for the mother spline
db0 = @(j,xi) motherSpline(j+4,xi);

switch dim,
  case 1, 
    TP(pad+(1:m)) = reshape(T,m,1);
    clear T;
    p = pad + p(1); xi = xi(1);
    
    Tc(valid) = TP(p+2).*b0(1,xi-2) + TP(p+1).*b0(2,xi-1) ...
      + TP(p).*b0(3,xi)+ TP(p-1).*b0(4,xi+1);
    if doDerivative,
      dT(valid) = TP(p+2).*db0(1,xi-2) + TP(p+1).*db0(2,xi-1) ...
        + TP(p).*db0(3,xi)+ TP(p-1).*db0(4,xi+1);
    end
  case 2, 
    TP(pad+(1:m(1)),pad+(1:m(2))) = T;
    clear T;
    p  = (pad + p(1)) + i2*(pad + p(2) - 1);
    for j1=-1:2,                              % Tc as weighted sum
      for j2=-1:2,
        Tc(valid) = Tc(valid) + TP(p+j1*i1+j2*i2).*b0(3-j1,xi(1)-j1).*b0(3-j2,xi(2)-j2);
        if doDerivative,
          dT(valid,1) = dT(valid,1) + ...
            TP(p+j1*i1+j2*i2).*db0(3-j1,xi(1)-j1).* b0(3-j2,xi(2)-j2);
          dT(valid,2) = dT(valid,2) + ...
            TP(p+j1*i1+j2*i2).* b0(3-j1,xi(1)-j1).*db0(3-j2,xi(2)-j2);
        end;
      end;
    end;
  case 3, 
    TP(pad+(1:m(1)),pad+(1:m(2)),pad+(1:m(3))) = T;
    clear T;
    p  = (pad + p(1)) + i2*(pad + p(2) - 1) + i3*(pad + p(3) -1);
    % compute Tc as weighted sum
    for j1=-1:2,
      for j2=-1:2,
        for j3=-1:2,
          Tc(valid) = Tc(valid) + TP(p+j1*i1+j2*i2+j3*i3).* ...
            b0(3-j1,xi(1)-j1).*b0(3-j2,xi(2)-j2).*b0(3-j3,xi(3)-j3);
          if doDerivative,
            dT(valid,1) = dT(valid,1) + TP(p+j1*i1+j2*i2+j3*i3).* ...
              db0(3-j1,xi(1)-j1).*b0(3-j2,xi(2)-j2).*b0(3-j3,xi(3)-j3);
            dT(valid,2) = dT(valid,2) + TP(p+j1*i1+j2*i2+j3*i3).* ...
              b0(3-j1,xi(1)-j1).*db0(3-j2,xi(2)-j2).*b0(3-j3,xi(3)-j3);
            dT(valid,3) = dT(valid,3) + TP(p+j1*i1+j2*i2+j3*i3).* ...
              b0(3-j1,xi(1)-j1).*b0(3-j2,xi(2)-j2).*db0(3-j3,xi(3)-j3);
          end;
        end;
      end;
    end;
end;

if doDerivative
    for i=1:dim, dT(:,i) = dT(:,i)/h(i); end
    if not(matrixFree)
        dT = spdiags(dT,n*(0:(dim-1)),n,dim*n);
    end
end
%==============================================================================

