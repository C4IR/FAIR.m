%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% conveniently sets the various parameters for the various optimization schemes
%
% schemes supported are
%     'SD'                steepestDescent
%     'Nesterov'          Nesterov
%     'Nesterov-GN'       NesterovGN
%     'GN'                GaussNewton, solver to be determined
%     'PIR-GN'            GaussNewton, solver = backslash
%     'NPIR-GN'           GaussNewton, solver = regularizer('get','solver')
%     'lBFGS'             lBFGS,        solver to be determined
%     'TrustRegion'       TrustRegion, solver to be determined
% =============================================================================

function opt = optPara(flag,varargin)

% administration, parameter setting, and such

if nargin == 0,
  help(mfilename);
  runMinimalExample;
  y = 'endOfMinimalExample';
  return;
end;

% default values for various schemes
opt.scheme      = '';
opt.maxIter     = 10;           % maximum number of iterations
opt.tolJ        = 1e-3;         % for stopping, objective function
opt.tolY        = 1e-2;         %   - " -     , current value
opt.tolG        = 1e-2;         %   - " -     , norm of gradient

switch flag,
  case {'SD','steepestDescent'},
    opt.scheme      = @steepestDescent;   % optimizer
    opt.lineSearch  = @Armijo;      % linesearch scheme
    opt.LSmaxIter   = 10;           % maximum number of line search iterations
    opt.LSreduction = 1e-4;         % minimal reduction in line search
    opt.stepLengt   = 1;

  case 'Nesterov',
    opt.scheme      = @Nesterov;    % optimizer
%    opt.lineSearch  = @Armijo;      % linesearch scheme
    opt.LSmaxIter   = 10;           % maximum number of line search iterations
    opt.LSreduction = 1e-4;         % minimal reduction in line search
    opt.lipschitzC   = 1;

  case 'Nesterov-GN',
    opt.scheme      = @NesterovGN;  % optimizer
    opt.solver      = 'backslash';  % solver for GN system
%    opt.lineSearch  = @Armijo;      % linesearch scheme
    opt.LSmaxIter   = 10;           % maximum number of line search iterations
    opt.LSreduction = 1e-4;         % minimal reduction in line search
    opt.lipschitzC   = 1;
    
  case 'GN',
    opt.scheme      = @GaussNewton; % optimizer
    opt.solver      = '';           % solver for linear system, see below
    opt.lineSearch  = @Armijo;      % linesearch scheme
    opt.LSmaxIter   = 10;           % maximum number of line search iterations
    opt.LSreduction = 1e-4;         % minimal reduction in line search

  case 'PIR-GN',
    opt.scheme      = @GaussNewton; % optimizer
    opt.solver      = 'backslash';  % solver for linear system, see below
    opt.lineSearch  = @Armijo;      % linesearch scheme
    opt.LSmaxIter   = 10;           % maximum number of line search iterations
    opt.LSreduction = 1e-4;         % minimal reduction in line search

  case 'NPIR-GN',
    opt.scheme      = @GaussNewton; % optimizer
    opt.solver      = '';           % solver for linear system, see below
    opt.lineSearch  = @Armijo;      % linesearch scheme
    opt.LSmaxIter   = 10;           % maximum number of line search iterations
    opt.LSreduction = 1e-4;         % minimal reduction in line search
    
    opt.solver = regularizer('get','solver');
    
  case 'lBFGS',
    opt.scheme      = @lBFGS;       % optimizer
    opt.solver      = '';           % solver for linear system, see below
    opt.lineSearch  = @Armijo;      % linesearch scheme
    opt.LSmaxIter   = 10;           % maximum number of line search iterations
    opt.LSreduction = 1e-4;         % minimal reduction in line search
    
    opt.solver = regularizer('get','solver');
    
  case 'TrustRegion',
    opt.scheme      = @TrustRegion; % optimizer
%    opt.solver      = '';           % solver for linear system, see below
%    opt.lineSearch  = @Armijo;      % linesearch scheme
%    opt.LSmaxIter   = 10;           % maximum number of line search iterations
%    opt.LSreduction = 1e-4;         % minimal reduction in line search
%    
    opt.solver = regularizer('get','solver');
    opt.preconditioner = @(x,para) x;
  otherwise
    keyboard
    error(flag)
end;

opt.Plots       = @FAIRplots;      % for plots;
% warning(sprintf('plots disabled in %s',mfilename))
% opt.Plots       = @(varargin) [];  % for no plots;
opt.vecNorm     = @norm;           % norm to be used for dJ and dy

% overwrite defaults
for k=1:2:length(varargin),     % overwrites default parameter
  eval(['opt.(''',varargin{k},''')=varargin{',int2str(k+1),'};']);
end;

%------------------------------------------------------------------------------
function runMinimalExample

E9_Hands_NPIR_OPT
%==============================================================================


