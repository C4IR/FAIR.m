%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% Main function for the image model, uses persistent parameter
% Typical call
%  1.   initialize (pick model and submit parameters)
%
%       imgModel('reset','imgModel','splineInter','regularizer','moments','theta',1e-2);
%
%  2.   use:        
%
%       [Tc,dT] = imgModel(coefT,omega,yc);
%      
%       evaluates the image model as a linear combination of basis functions with 
%       the coeeficients coeffT of size size(coefT) at a cell-centered grid of points
%       y = [Y(:,1);Y(:,2),...,Y(:,d)] using parameters specified by the persistent variable OPTN
% 
% 
%       [Tc,dT] = imgModel(dataT,omega,yc,specific{:});
%       overwrites the persistent options (not permanently)
%
%  3.   additional option 'coefficients':
%
%       [coefT,coefR] = imgModel('coefficients',dataT,dataR,omega);
%        
%       returns the coefficients of a basis representation, note
%       coefT = dataT                     for linear interpolation
%       coefT = getSplineCoefficients     for spline interpolation
% =======================================================================================

function varargout = imgModel(varargin)

% handle options
persistent OPTN 

if nargin == 0 & nargout == 0 & isempty(OPTN),
  help(mfilename);
  runMinimalExample;
  return;
end;

% check for reset, set, clear, disp, see dealOptions.m
[method,OPTN,task,stop] = dealOptions(OPTN,varargin{:});
if stop,  varargout = {method,OPTN};  return; end

% update spline coefficients, if neccessary
if strcmp(task,'coefficients'),
  T   = varargin{2};  
  R   = varargin{3};  
  dim = length(varargin{4})/2;
  if ~isempty(strfind(method,'spline')),
    T = getSplineCoefficients(T,'dim',dim,OPTN{:},varargin{5:end});
    R = getSplineCoefficients(R,'dim',dim,OPTN{:},varargin{5:end});
  end;
  varargout = {T,R};
  return;
end;

% do the work
[method,optn] = dealOptions(OPTN,'set','doDerivative',(nargout>1),varargin{4:end});
k = find(strcmp(optn,mfilename));
optn (k:k+1) = [];
T         = varargin{1};
omega     = varargin{2};
x         = varargin{3}(:); 
[T,dT]    = feval(method,T,omega,x,optn{:});
varargout = {T,dT};

%------------------------------------------------------------------------------

function runMinimalExample
fprintf('run minimal example for %s\n',mfilename)
% load data
setup2DhandData
%==============================================================================

