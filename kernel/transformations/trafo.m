%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function varargout = trafo(varargin)
%
% Main function for the transformation model, uses persistent parameter
% Typical call
%    [y,dy] = trafo(w,X);
%    computes y = Q(X)*f(w), dy = Q(X)df(w), 
%    where Q(X) is a sample of the basis functions on the grid X,
%       Q is allocated persistent in the method to be used, X is used for proper size
%
%  1.     initialize (pick model and submit parameters)
%
%         trafo('reset','trafo','rigid2D','center',center);
%
%  2.     [y,dy] = trafo(w,X);
%        evaluates the transformation model, e.g. yc = Q(x)*f(w) and its derivative
%
%  3.   w0 = trafo('w0'); 
%        returns parameterization of identity
% 
%  4.   [y,dy] = trafo(w,x,specific{:}), uses specific parameters in this call
%
% see also E4_US_trafo, transformations/contents.m.
%
% =======================================================================================

function varargout = trafo(varargin)

persistent OPTN 

if nargin == 0 && nargout == 0 && isempty(OPTN),
  help(mfilename);
  runMinimalExample;
  return;
end;

% handle options
[method,OPTN,task,stop] = dealOptions(OPTN,varargin{:});
if stop,
  varargout{1} = method;
  if any(strcmp(task,{'reset','set'})), trafo('w0'); end; % clear Q
  if nargout > 1, varargout{2} = OPTN;  end;
  return;
end

% return parameters for identity transformation
if strcmp(task, 'w0'),
  [~,w0] = feval(method,'w0',[],OPTN{:});
  varargout = {w0};
  return;
end;

% do the work
[method,optn] = dealOptions(OPTN,'set','doDerivative',(nargout>1),varargin{3:end});
doDerivative  = dealOptions(optn,'get','doDerivative');
w = varargin{1};
x = varargin{2};
if ~doDerivative,
  y  = feval(method,w,x,optn{:});
  varargout = {y,[]};
  return;
end;


[y,dy] = feval(method,w,x,optn{:});
if ~strcmp(dealOptions(optn,'get','debug'),'on'), 
  varargout = {y,dy};
  return; 
end;

% DEBUGGING mode (used for derivative checks),
% transform sparse presentation of dy to full
if iscell(dy)  
  if length(dy) == 1,
    p = size(dy{1});
    dim = numel(w)/p(2);
    dy = kron(speye(dim),dy{1});
  else
    error('nyi')
  end;
end;
varargout = {y,dy};

%------------------------------------------------------------------------------

function runMinimalExample
fprintf('run minimal example for %s\n',mfilename)
E4_US_rotation
%==============================================================================
