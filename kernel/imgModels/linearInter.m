%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
% ##2
%==============================================================================
% function [Tc,dT] = linearInter(T,omega,x,varargin)
%
% linear interpolator for the data T given on a cell-centered grid 
% to be evaluated at locations x
%
% Input
%  T        coefficients for image model img : \R^{dim*n} \to\R^{n)}, 
%           for j=1:n, img(j)  = sum T_i b_i(x_j); end;
%  omega    specification of domain
%           Omega = (omega(1),omega(2)) x ... x (omega(2*dim-1,omega(2*dim)), 
%           spatial dimension dim = length(omega/2)
%  x        location evaluation points
%           x = reshape(x,n,dim), x(j,:) coordinates of j-th point
%  varargin additional parameters such as doDerivative,         
%##4 cleanup        
%==============================================================================

function [Tc,dT] = linearInter(T,omega,x,varargin)
         
Tc = mfilename('fullpath'); dT = []; 

if nargin == 0, 
  help(mfilename);
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

%------------------------------------------------------------------------------

% get data size m, cell size h, dimension d, and number n of interpolation points
dim = length(omega)/2;
m   = size(T);         if dim == 1, m = numel(T); end;
h   = (omega(2:2:end)-omega(1:2:end))./m;
n   = length(x)/dim;
x   = reshape(x,n,dim);

% map x from [h/2,omega-h/2] -> [1,m],
for i=1:dim, x(:,i) = (x(:,i)-omega(2*i-1))/h(i) + 0.5; end;

Tc = zeros(n,1); dT = [];                     % initialize output
if doDerivative, dT = zeros(n,dim);  end;     % allocate memory in column format
Valid = @(j) (0<x(:,j) & x(:,j)<m(j)+1);      % determine indices of valid points

switch dim,
  case 1, valid = find( Valid(1) );   
  case 2, valid = find( Valid(1) & Valid(2) );   
  case 3, valid = find( Valid(1) & Valid(2) & Valid(3) );   
end;

if isempty(valid),                        
  if doDerivative, dT = sparse(n,dim*n); end; % allocate memory in column format
  return; 
end;

pad = 1; TP = zeros(m+2*pad);                 % pad data to reduce cases

P = floor(x); x = x-P;                        % split x into integer/remainder
p = @(j) P(valid,j); xi = @(j) x(valid,j);

% increments for linearized ordering
i1 = 1; i2 = size(T,1)+2*pad; i3 = (size(T,1)+2*pad)*(size(T,2)+2*pad);

switch dim,
  case 1, 
    TP(pad+(1:m)) = reshape(T,m,1);
    clear T;
    p = pad + p(1);
    Tc(valid) = TP(p).* (1-xi(1)) + TP(p+1).*xi(1);   % compute weighted sum
    
    if ~doDerivative, return; end;             
    % compute and format the derivative
    dT(valid) = TP(p+1)-TP(p);
  case 2, 
    TP(pad+(1:m(1)),pad+(1:m(2))) = T;
    clear T;
    p  = (pad + p(1)) + i2*(pad + p(2) - 1);
    % compute Tc as weighted sum
    Tc(valid) = (TP(p)   .* (1-xi(1)) + TP(p+i1)    .*xi(1)) .* (1-xi(2)) ...
      + (TP(p+i2) .* (1-xi(1)) + TP(p+i1+i2) .*xi(1)) .* (xi(2));
  
    if ~doDerivative, return; end;
    dT(valid,1) = (TP(p+i1)-TP(p)).*(1-xi(2)) + (TP(p+i1+i2)-TP(p+i2)).*xi(2);
    dT(valid,2) = (TP(p+i2)-TP(p)).*(1-xi(1)) + (TP(p+i1+i2)-TP(p+i1)).*xi(1);
  case 3, 
    TP(pad+(1:m(1)),pad+(1:m(2)),pad+(1:m(3))) = T;
    clear T;
    p  = (pad + p(1)) + i2*(pad + p(2) - 1) + i3*(pad + p(3) -1);
    % compute Tc as weighted sum
    Tc(valid) = ((TP(p).*(1-xi(1))+TP(p+i1).*xi(1)).*(1-xi(2))...
      +(TP(p+i2).*(1-xi(1))+TP(p+i1+i2).*xi(1)).*(xi(2))).*(1-xi(3)) ...
      +((TP(p+i3).*(1-xi(1))+TP(p+i1+i3).*xi(1)).*(1-xi(2)) ...
      +(TP(p+i2+i3).*(1-xi(1))+TP(p+i1+i2+i3).*xi(1)).*(xi(2))).*(xi(3));
    
    if ~doDerivative, return; end;
    dT(valid,1) = ((TP(p+i1)-TP(p)).*(1-xi(2))+(TP(p+i1+i2)-TP(p+i2)).*xi(2)).*(1-xi(3)) ...
      +((TP(p+i1+i3)-TP(p+i3)).*(1-xi(2))+(TP(p+i1+i2+i3)-TP(p+i2+i3)).*xi(2)).*(xi(3));
    dT(valid,2) = ((TP(p+i2)-TP(p)).*(1-xi(1))+(TP(p+i1+i2)-TP(p+i1)).*xi(1)).*(1-xi(3)) ...
      +((TP(p+i2+i3)-TP(p+i3)).*(1-xi(1))+(TP(p+i1+i2+i3)-TP(p+i1+i3)).*xi(1)).*(xi(3));
    dT(valid,3) = ((TP(p+i3).*(1-xi(1))+TP(p+i1+i3).*xi(1)).*(1-xi(2)) ...
      +(TP(p+i2+i3).*(1-xi(1))+TP(p+i1+i2+i3).*xi(1)).*(xi(2))) ....
      -((TP(p).*(1-xi(1))+TP(p+i1).*xi(1)).*(1-xi(2)) ...
      +(TP(p+i2).*(1-xi(1))+TP(p+i1+i2).*xi(1)).*(xi(2)));
end;
if doDerivative
    for i=1:dim, dT(:,i) = dT(:,i)/h(i); end
    if not(matrixFree)
        dT = spdiags(dT,n*(0:(dim-1)),n,dim*n);
    end
end
%==============================================================================
