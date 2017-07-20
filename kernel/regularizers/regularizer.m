%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function varargout = regularizer(varargin)
%
% Main function for the regularization model, uses persistent parameter
%
% The regularization functional is either based on a linear differential
% operator B, i.e.
%
%        S(Y) = 0.5*alpha*hd*|B*Y|^2
%
%  or a general (non-linear) functional S(Y)
%
% Typical call
%    [Sc,dS,d2S] = regularizer(Y,omega,m)
%
%  1.     initialize (pick regularizer and set parameters)
%
%         regularizer('reset','regularizer','mfElastic','alpha',1e3,...
%               'mu',1,'lambda',0','rigid2D','center',center);
%
%  2.     [Sc,dS,d2S] = regularizer(Y,omega,m)
%        evaluates the regularization model, e.g. Sc = 0.5*yc'*B'*B*yc and its derivatives
%
%  A = alpha*hd*B'*B is made persistent for eficiency
%  hd = prod(omega./m)
%  an option grid is introduce to indicate the appropriate discretization:
%  staggered for elastic and diffusive,
%  cell-centered for curvature
%  nodal for TV, hyper elastic or diffusive EPI regularization
%==============================================================================

function varargout = regularizer(varargin)

persistent OPTN A

if nargin == 0 && nargout == 0 && isempty(OPTN),
    help(mfilename);
    return;
end;

% -----------------------------------------------------------------------------
% handle options
[method,OPTN,task,stop] = dealOptions(OPTN,varargin{:});


% check the grid according to the regularizer to be used
% setup default solver for Gauss-Newton systems

if strcmp(task,'set') || strcmp(task,'reset'),
    switch method,
        case 'mfElastic',
            scheme     = 'elastic';
            matrixFree = 1;
            grid       = 'staggered';
            solver     = 'MG-elastic';
            
        case 'mbElastic',
            scheme     = 'elastic';
            matrixFree = 0;
            grid       = 'staggered';
            solver     = 'backslash';
            
        case 'mbElasticNodal',
            scheme     = 'elasticNodal';
            matrixFree = 0;
            grid       = 'nodal';
            solver     = 'backslash';
            
        case 'mbCurvature',
            scheme     = 'curvature';
            matrixFree = 0;
            grid       = 'cell-centered';
            solver     = 'backslash';
            
        case 'mfCurvature',
            scheme     = 'curvature';
            matrixFree = 1;
            grid       = 'cell-centered';
            solver     = 'PCG-curvature';
            
        case 'mbHyperElastic',
            scheme     = 'hyperElastic';
            matrixFree = 0;
            grid       = 'nodal';
            solver     = 'backslash';
            
            
        case 'mfHyperElastic',
            scheme     = 'hyperElastic';
            matrixFree = 1;
            grid       = 'nodal';
            solver     = 'PCG-hyperElastic';
            
        case 'mbHyperElasticFEM',
            scheme     = 'hyperElasticFEM';
            matrixFree = 0;
            grid       = 'FEM';
            solver     = 'backslash';
            
        case 'mfHyperElasticFEM',
            scheme     = 'hyperElasticFEM';
            matrixFree = 1;
            grid       = 'FEM';
            solver     = 'PCG-hyperElastic';

        case 'mbElasticFEM',
            scheme     = 'elasticFEM';
            matrixFree = 0;
            grid       = 'FEM';
            solver     = 'backslash';
            
        otherwise
            scheme = method;
            [grid,matrixFree,solver] = feval(scheme,'para',[],[],OPTN{:});
    end;
    
    [dummy,OPTN] = dealOptions(OPTN,'set','scheme',scheme,...
        'grid',grid,'matrixFree',matrixFree,'solver',solver);
end;

% return, if no further work has to be handled
if stop,
    varargout{1} = method;
    if nargout > 1, varargout{2} = OPTN;  varargout{3} = A;  end;
    return;
end
% -----------------------------------------------------------------------------
% do the work

% extract regularization parameters
scheme      = dealOptions(OPTN,'get','scheme');
alpha       = dealOptions(OPTN,'get','alpha');
matrixFree  = dealOptions(OPTN,'get','matrixFree');

% extract variables
yc     = varargin{1};
omega  = varargin{2};
m      = varargin{3};
doDerivative = (nargout>1);
varargin = varargin(4:end);

for k=1:2:length(varargin), % overwrites defaults
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;
[Sc,dS,d2S] = feval(scheme,yc,omega,m,...
    'alpha',alpha,'matrixFree',matrixFree,...
    'doDerivative',doDerivative,OPTN{:},varargin{:});
varargout = {Sc,dS,d2S};

%------------------------------------------------------------------------------

function runMinimalExample
fprintf('run minimal example for %s\n',mfilename)
keyboard
%==============================================================================


