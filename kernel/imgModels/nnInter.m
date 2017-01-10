%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
% ##2
%==============================================================================
%
% function [Tc,dT] = nnInter(T,omega,x,varargin);
%
% next neighbour interpolator for the data T given on a cell-cantered grid 
% evaluated at x, see FAIR 3.2 p24.
%==============================================================================

function [Tc,dT] = nnInter(T,omega,x,varargin)
         
Tc = mfilename('fullpath'); dT = []; 

if nargin == 0, 
    
  testOneImgModel(mfilename);
  return;
elseif nargin == 1 && isempty(T),
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
n   = length(x)/dim;  
dT  = sparse(n,dim*n);
x   = reshape(x,n,dim);
% map x from [h/2,omega-h/2] -> [1,m],
for i=1:dim, x(:,i) = (x(:,i)-omega(2*i-1))/h(i) + 0.5; end;

% round x
x = round(x);

Tc = zeros(n,1);                            % initialize output
Valid = @(j) (0<x(:,j) & x(:,j)<m(j)+1);    % determine indices of valid points
switch dim,
  case 1, valid = find( Valid(1) );   
  case 2, valid = find( Valid(1) & Valid(2) );   
  case 3, valid = find( Valid(1) & Valid(2) & Valid(3) );   
end;
             
if isempty(valid),                        
  return; 
end;

p = @(j) x(valid,j);
% increments for linearized ordering
i1 = 1; i2 = size(T,1); i3 = size(T,1)*size(T,2);
switch dim,
  case 1, p = p(1);
  case 2, p = p(1) + i2*(p(2) - 1);
  case 3, p = p(1) + i2*(p(2) - 1) + i3*(p(3) -1);
end;
Tc(valid) = T(p);
%==============================================================================
