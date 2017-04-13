%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% function [fh,ph,th] = checkDerivative(fctn,x0);
%
% checks the implementation of a derivative by comparing the function 
% with the Taylor-poly
%
%   \| f(x0 + h ) - TP_p(x0,h) \|   !=   O( h^{p+1} ) 
%
% Input:
%  fctn    function handle
%  x0      expanding point
%
% Output:
%  fh      figure handle to graphical output
%  ph      plot handle
%  th      text handle
% call checkDerivative for a minimal example.
%==============================================================================

function varargout = checkDerivative(fctn,x0,varargin)

if nargin == 0, % help and minimal example
  help(mfilename); 
  fctn = @xSquare;  
  x0   = 1;  
  checkDerivative(fctn,x0);
  return;
end;


fig = [];
for k=1:2:length(varargin),     % overwrites default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

FAIRmessage(mfilename,'.')
fprintf('test derivative of function <%s>\n',func2str(fctn));

fprintf('T0 = |f0 - ft|, T1 = |f0+h*f0'' - ft|\n');
[f0,df] = feval(fctn,x0);  
if ~isnumeric(df),
  try
    [f0,para,df] = feval(fctn,x0);  
  catch
    if isfield(df,'Q'), df = df.Q;  end;
  end;
end;

h = logspace(-1,-10,10);
v = randn(size(x0)); 
if isnumeric(df),
  dvf = df*v; 
elseif isa(df,'function_handle'),
  dvf = df(v);
else
  keyboard;
end;

for j=1:length(h),
  ft = feval(fctn,x0+h(j)*v);      % function value
  T0(j) = norm(f0-ft);             % TaylorPoly 0
  T1(j) = norm(f0+h(j)*dvf - ft);  % TaylorPoly 1
  fprintf('h=%12.4e     T0=%12.4e    T1=%12.4e\n',h(j),T0(j),T1(j));
end;

fh = FAIRfigure(fig);
ph = loglog(h,[T0;T1]); set(ph(2),'linestyle','--')
th = title(sprintf('%s: |f-f(h)|,|f+h*dvf -f(h)| vs. h',mfilename));

if nargout>0,
  varargout = {fh,ph,th};
end;
FAIRmessage('.');
%------------------------------------------------------------------------------
function [y,dy] = xSquare(x)
if nargin == 0, return; end;
y = x.^2; dy = reshape(2*x,1,[]);
%==============================================================================
