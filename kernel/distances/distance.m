%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function varargout = distance(varargin)
%
% Main function for distance measures, possible calls:
%
% initialize: distance('reset','distance','SSD');
% use:        [Dc,rc,dD,drc,d2psi] = distance(Tc,Rc,omega,m);
% evaluates the distances between Tc and Rc, cell volume (pixel/voxel) 
% is based on omega and m;
% internals: residual rc and outer function psi, D(Tc) = psi(rc(Tc))
%
% administration, parameter setting:
%  for resetting, intitializing, setting, updating, clearing, displaying,
%  options - deals parameterization
%
% specific options:
%  [Dc,rc,dD,drc,d2psi] = distance(Tc,Rc,omega,m,specific{:});
%
% see also E9_Hands_MLIR_SSD_mbElas
% =======================================================================================

function varargout = distance(varargin);

persistent OPTN

if nargin == 0 && nargout == 0 && isempty(OPTN),
  help(mfilename);
  return;
end;

% handle options
[method,OPTN,task,stop] = dealOptions(OPTN,varargin{:});
if stop,
  varargout{1} = method;
  if nargout > 1, varargout{2} = OPTN;  end;
  return;
end

% do the work

Tc      = varargin{1};
Rc      = varargin{2};
omega   = varargin{3};
m       = varargin{4};

% check if weights are to be applied and extract current level
k = find(strcmp(OPTN(1:2:end),'weights'));
if not(isempty(k))
    MLw = OPTN{2*k};
    for k=1:length(MLw)
        if isfield(MLw{k},'m') && all(MLw{k}.m == m)
            Wc = MLw{k}.Wc;
            break
        end;
    end;
    varargin{end+1} = 'weights';
    varargin{end+1} = Wc;
end;

[method,optn]  = dealOptions(OPTN,'set','doDerivative',(nargout>3),varargin{5:end});
[Dc,rc,dD,dr,d2psi] =  feval(method,Tc,Rc,omega,m,optn{:});
  
varargout = {Dc,rc,dD,dr,d2psi};
%==============================================================================

