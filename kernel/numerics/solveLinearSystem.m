%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration.
% For details see
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% This function manages the solution of H * dy = - dJ' for various settings
%
%   'backslash',         simply solve dy = H\rhs;
%   'mbCG',               matrixBased, plain CG
%   'mbPCG-SGS',         matrixBased, PCG with Symmetric Gauss Seidel preconditioner
%   'mbPCG-ICHOL',     matrixBased, PCG with incomplete Cholesky preconditioner
%   'mbPCG-Jacobi', matrixBased, PCG with Jacobi preconditioner
%
%   'MG-elastic',           matrixFree, PCG Multigrid for elastic regularizer
%   'CG-elastic',       plain CG with H = A as in (1) below
%   'CG-curvature',          plain CG with H = A as in (1) below
%   'PCG-elastic',      Jacobi preconditioned with H = A as in (1,) below
%   'PCG-curvature',      Jacobi preconditioned with H = A as in (1,) below
%   'PCG-hyperElastic', ???
%   'CG'                      plain CG with H full or as operator
%
% the operator H can be a
%    - matrix (matrixBased)
%    - function (coding the action ofH)
%    - struct where H describes the pieces of a complex H
%
%  (1) action of H: distance + regularizer,
%       H = P'*dr'*d2psi*dr*P + d2S
%
%      Hoperator = @(x) ...
%       H.d2D.P((H.d2D.dr'*H.d2D.d2psi*H.d2D.dr)*H.d2D.P(x)) ...
%       + H.d2S.d2S(x,H.omega,H.m);
%
%  (2) Jacobi (diagonal) preconditing
%      Ddiag = diag(H.d2D.dr'*H.d2D.d2psi*H.d2D.dr);
%      D = H.d2D.P(full(Ddiag))  +  H.d2S.diag(H.omega,H.m);
%      PC = @(x) D.\x; % Jacobi preconditioner
%==============================================================================

function [dy,solver] = solveLinearSystem(rhs,H,solver,varargin)

if nargin == 0,
  testSolver
  help(mfilename);
  runMinimalExample;
  dy = 'endOfMinimalExample';
  return;
end;

maxIterCG = 500;
tolCG     = 1e-1;

for k=1:2:length(varargin), % overwrites default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

if isempty(solver) && isnumeric(H),
  warning('no solver specified, try backslash');
  solver = 'backslash';
end;

if isa(solver,'function_handle')
  dy     = solver(rhs,H,maxIterCG,tolCG);
  return;
end

if isstruct(H), % matrixFree mode, configure operator
  Hoperator = @(x) ...
    H.d2D.P((H.d2D.dr'*H.d2D.d2psi*H.d2D.dr)*H.d2D.P(x)) ...
    + H.d2S.d2S(x,H.omega,H.m);
end;

%------------------------------------------------------------------------------
switch solver
  %------------------------------------------------------------------------------
  
  % ---------------------------------------------------------------------------
  % matrix based
  % ---------------------------------------------------------------------------
  
  case 'backslash',
    dy = H\rhs;
    
  case 'mbCG',
    [dy,flag,relres,iter] = pcg(H,rhs,tolCG,maxIterCG);
    
  case 'mbPCG-SGS',
    L   = tril(H); % Symmetric Gauss Seidel Preconditioning,
    D   = diag(H); % L is lower, D is diagonal, U = L'
    SGS = @(x) L\(D.*(L'\x));
    [dy,flag,relres,iter] = pcg(H,rhs,tolCG,maxIterCG,SGS);
    
  case 'mbPCG-ICHOL',
    %L1   = cholinc(sparse(H),'0');
    L = ichol(H);
    [dy,flag,relres,iter] = pcg(H,rhs,tolCG,maxIterCG,L,L');
    
  case 'mbPCG-Jacobi',
    D   = diag(H); % D is diagonal
    PC = @(x) D.\x;
    [dy,flag,relres,iter] = pcg(H,rhs,tolCG,maxIterCG,PC);
    
    
    
    % ---------------------------------------------------------------------------
    % matrix free
    % ---------------------------------------------------------------------------
    
    % In all matrix free solvers the approximates Hessian is supplied by function handle.
    % Note, that the functionals have the common form
    %
    %   Jc = D + S (+ P) ==> d_2 Jc = d2D + d2S (+d2P)
    %
    % regularization and penalization (S,P) supply a struct stored in H.d2S or H.d2P
    % which already has a function handle describing the operator.
    % The  distance term (d2D) needs more care and the approximation to the Hessian
    % depends on the solver as well. Therefore the objective function does NOT supply
    % the function handle, but all variables needed to built such.
    % For multi-grid (available for elastic and curvature), only the diagonal of the distance
    % term is used
    % whereas in conjugate gradient (CG, available for elastic, curvature and hyperelastic)
    % methods we can model also the off-diagonals d2D(x) = (P'*dr'*d2psi*dr*P) * x
    %
    % The default solver for a chosen regularizer is parameterized by d2S.solver.
    % However, one can overload this setting using
    %    >> MLIR(..., 'solverNPIR','myFavoriteSolver', ... );
    
  case {'MG-elastic','MG'}, % multigrid for elastic regularizer
    dy = MGsolver(rhs,H);
    
  case {'CG-elastic','CG-curvature'},
    if isstruct(H)
      %       A         = @(x) ...
      %         H.d2D.P((H.d2D.dr'*H.d2D.d2psi*H.d2D.dr)*H.d2D.P(x)) ...
      %         + H.d2S.d2S(x,H.omega,H.m);
      [dy,flag,relres,iter] = pcg(Hoperator,rhs,tolCG,maxIterCG);
    else
      [dy,flag,relres,iter] = pcg(H,rhs,tolCG,maxIterCG);
    end
    
  case {'PCG-elastic','PCG-curvature'}
    
    % operator
    %     A         = @(x) ...
    %       H.d2D.P((H.d2D.dr'*H.d2D.d2psi*H.d2D.dr)*H.d2D.P(x)) ...
    %       + H.d2S.d2S(x,H.omega,H.m);
    % preconditioner
    Ddiag = diag(H.d2D.dr'*H.d2D.d2psi*H.d2D.dr);
    D = H.d2D.P(full(Ddiag))  +  H.d2S.diag(H.omega,H.m);
    PC = @(x) D.\x; % Jacobi preconditioner
    [dy,flag,relres,iter] = pcg(Hoperator,rhs,tolCG,maxIterCG,PC);
    
  case 'CG'
    if isstruct(H)
      [dy,flag,relres,iter] = pcg(H.operator,rhs,tolCG,maxIterCG);
    else
      [dy,flag,relres,iter] = pcg(H,rhs,tolCG,maxIterCG);
    end
    
  case {'PCG-hyperElastic'}
    M         = @(x) H.d2D.P((H.d2D.dr'*H.d2D.d2psi*H.d2D.dr)*H.d2D.P(x));
    Hoperator = @(x) M(x) + H.d2S.d2S(x,H.omega,H.m,H.d2S.yc);
    Ddiag     = diag(H.d2D.dr'*H.d2D.d2psi*H.d2D.dr);
    D         = H.d2D.P(full(Ddiag))  +  H.d2S.diag(H.d2S.yc);
    Preconditioner = @(x) D.\x; % Jacobi preconditioner
    [dy,flag,relres,iter] = pcg(Hoperator,rhs,tolCG,maxIterCG,Preconditioner);
    
  otherwise,
    keyboard
    error(1)
    
    %------------------------------------------------------------------------------
end
%------------------------------------------------------------------------------

if exist('flag','var')
  flags = {
    sprintf('iter=%d>maxIter=%d but relres=%e>relResTol=%s',...
    iter,maxIterCG,relres,tolCG)
    'preconditioner is ill-conditioned'
    'stagnation: x_k=x_{k+1}'
    'scalars too small/large to continue computation'
    };
  switch flag
    case {1,2,3,4}
      fprintf('%s // PCG: %s\n',flags{flag});
    otherwise
      fprintf('%s // PCG :iter = %d of %d, relRes= %1.2e, tolCG= %1.2e\n',...
        mfilename,iter,maxIterCG,relres,tolCG);
  end
end

%------------------------------------------------------------------------------

function runMinimalExample


%% prepare 3D data
setup3DbrainData
level = 4;
omega = ML{level}.omega;
m     = ML{level}.m;

viewImage('reset','viewImage','viewImage2D','colormap',bone(256),'axis','off');
imgModel('reset','imgModel','splineInterMex','regularizer','moments','theta',1e-2);
distance('reset','distance','SSD');
trafo('reset','trafo','affine3Dsparse');

[T,R] = imgModel('coefficients',ML{level}.T,ML{level}.R,omega,'out',0);
xc = getCellCenteredGrid(omega,m);
Rc = imgModel(R,omega,xc);
w0   = trafo('w0');
beta = 0; M = []; wRef = [];
[Jc,para,dJPIRsparse,HPIRsparse] = PIRobjFctn(T,Rc,omega,m,beta,M,wRef,xc,w0);


%% prepare 2D data and approximations to Hessian
setup2DhandData

viewImage('reset','viewImage','viewImage2D','colormap',bone(256),'axis','off');
imgModel('reset','imgModel','splineInter','regularizer','moments','theta',1e-2);
distance('reset','distance','SSD');

level = 4;
omega = ML{level}.omega;
m     = ML{level}.m;
[T,R] = imgModel('coefficients',ML{level}.T,ML{level}.R,omega,'out',0);

% initialize the image model, distance measure and regularizer
xc    = getCellCenteredGrid(omega,m);
Rc    = imgModel(R,omega,xc);

% setup PIR objective
w0   = trafo('w0');
beta = 1e3;
M    = diag([4,4,4,4,1,1]);

trafo('reset','trafo','affine2D');
w0   = trafo('w0');
wRef = w0;
[Jc,para,dJPIR,HPIR] = PIRobjFctn(T,Rc,omega,m,beta,M,wRef,xc,w0);

% setup NPIR objective
yc   = getStaggeredGrid(omega,m);
yRef = 0*yc;

regularizer('reset','regularizer','mbElastic','alpha',1e4,'mu',1,'lambda',0);
[Jc,para,dJmbElas,HmbElas] = NPIRobjFctn(T,Rc,omega,m,yRef,yc);

regularizer('reset','regularizer','mfElastic','alpha',1e4,'mu',1,'lambda',0);
[Jc,para,dJmfElas,HmfElas] = NPIRobjFctn(T,Rc,omega,m,yRef,yc);

HxElas = @(x) HmfElas.d2D.P((HmfElas.d2D.dr'*HmfElas.d2D.d2psi*HmfElas.d2D.dr)*HmfElas.d2D.P(x)) ...
  + HmfElas.d2S.d2S(x,HmfElas.omega,HmfElas.m);

yc   = getCellCenteredGrid(omega,m);
yRef = 0*yc;

regularizer('reset','regularizer','mbCurvature','alpha',1e4);
[Jc,para,dJmbCurv,HmbCurv] = NPIRobjFctn(T,Rc,omega,m,yRef,yc);

regularizer('reset','regularizer','mfCurvature','alpha',1e4);
[Jc,para,dJmfCurv,HmfCurv] = NPIRobjFctn(T,Rc,omega,m,yRef,yc);
HxCurv = @(x) HmfCurv.d2D.P((HmfCurv.d2D.dr'*HmfCurv.d2D.d2psi*HmfCurv.d2D.dr)*HmfCurv.d2D.P(x)) ...
  + HmfCurv.d2S.d2S(x,HmfCurv.omega,HmfCurv.m);

problems = {
  'PIR-backslash'
  'PIR-sparse'
  'mbElastic-backslash'
  'mbElastic-PCG-SGS'
  'mbElastic-PCG-ICHOL'
  'mbElastic-PCG-Jacobi'
  'mbElastic-CG'
  'mfElastic-MG'
  'mfElastic-CG'
  'mfElastic-PCG'
  'mbCurvature-backslash'
  'mbCurvature-PCG-SGS'
  'mbCurvature-PCG-ICHOL'
  'mbCurvature-PCG-Jacobi'
  'mbCurvature-CG'
  'mfCurvature-PCG'
  'mfCurvature-CG'
  }

for j = 1:length(problems)
  
  para = {};
  
  switch problems{j},
    
    case 'PIR-backslash',
      solver = 'backslash';
      H = HPIR; b = -dJPIR'; Hx = @(y) H*y;
      
    case 'PIR-sparse',
      solver = 'backslash';
      H = HPIRsparse; b = -dJPIRsparse'; Hx = @(y) H*y;
      
      %     case 2, % PRIR with conensed matrix
      % % fprintf('test PIR with condensed matrix\n');
      % % trafo('reset','trafo','affine2Dsparse');
      % % [Jc,para,dJ,H] = PIRobjFctn(T,Rc,omega,m,beta,M,wRef,xc,w0)
      % % dy = solveLinearSystem(-dJ',H,'backslash');
      % % test2 = norm(H*dy + dJ')
      
    case 'mbElastic-backslash'
      solver = 'backslash';
      H = HmbElas; b = -dJmbElas'; Hx = @(y) H*y;
      
    case 'mbElastic-PCG-SGS',
      solver = 'mbPCG-SGS';
      H = HmbElas; b = -dJmbElas'; Hx = @(y) H*y;
      para = {'tolCG',1e-5};
      
    case 'mbElastic-PCG-ICHOL',
      solver = 'mbPCG-ICHOL';
      H = HmbElas; b = -dJmbElas'; Hx = @(y) H*y;
      para = {'tolCG',1e-5};
      
    case 'mbElastic-PCG-Jacobi'
      solver = 'mbPCG-Jacobi';
      H = HmbElas; b = -dJmbElas'; Hx = @(y) H*y;
      para = {'tolCG',1e-5};
      
    case 'mbElastic-CG'
      solver = 'mbCG';
      H = HmbElas; b = -dJmbElas'; Hx = @(y) H*y;
      para = {'tolCG',1e-5};
      
    case 'mfElastic-MG'
      solver = 'MG-elastic';
      H = HmfElas; b = -dJmfElas'; Hx = HxElas;
      
    case 'mfElastic-CG'
      solver = 'CG-elastic';
      H = HmfElas; b = -dJmfElas'; Hx = HxElas;
      
    case 'mfElastic-PCG'
      solver = 'PCG-elastic';
      H = HmfElas; b = -dJmfElas'; Hx = HxElas;
      
    case 'mbCurvature-backslash'
      solver = 'backslash';
      H = HmbElas; b = -dJmbElas'; Hx = @(y) H*y;
      
    case 'mbCurvature-PCG-SGS',
      solver = 'mbPCG-SGS';
      H = HmbElas; b = -dJmbElas'; Hx = @(y) H*y;
      para = {'tolCG',1e-5};
      
    case 'mbCurvature-PCG-ICHOL',
      solver = 'mbPCG-ICHOL';
      H = HmbElas; b = -dJmbElas'; Hx = @(y) H*y;
      para = {'tolCG',1e-5};
      
    case 'mbCurvature-PCG-Jacobi'
      solver = 'mbPCG-Jacobi';
      H = HmbElas; b = -dJmbElas'; Hx = @(y) H*y;
      para = {'tolCG',1e-5};
      
    case 'mbCurvature-CG'
      solver = 'mbCG';
      H = HmbElas; b = -dJmbElas'; Hx = @(y) H*y;
      para = {'tolCG',1e-5};
      
      
    case 'mfCurvature-PCG'
      solver = 'PCG-curvature';
      H = HmfElas; b = -dJmfElas'; Hx = HxElas;
      
    case 'mfCurvature-CG'
      solver = 'CG-curvature';
      H = HmfElas; b = -dJmfElas'; Hx = HxElas;
      
    case 'MG'
      
    otherwise,
      error('1');
  end;
  
  dy = solveLinearSystem(b,H,solver,para{:});
  test = norm(Hx(dy)-b)/norm(b);
  fprintf('%-2d of %2d: %-40s test = %s\n',j,length(problems),...
    sprintf('test <%s>:',problems{j}),num2str(test));
end;

return;
%==============================================================================

%   case 'mfCurvature',
%     Afun     = @(dy) mfAy(dy,H);
%     [dy,FLAG] = pcg(Afun,rhs,tolCG,maxIterCG);
%     %

%
%   if isnumeric(H),
%
%   elseif isstruct(H),
%     if isfield(H,'solver'),
%       solver = H.solver;
%     elseif isfield(H,'d2S') && isfield(H.d2S,'solver'),
%       solver = H.d2S.solver;
%     else
%       error('solver has not been defined')
%     end;
%   else
%     keyboard
%   end;
%   fprintf('[set solver to <%s> in %s]\n',solver,mfilename);
% end;
% if isempty(solver)
%   if isstruct(H),
%   else % no regularizer initialized, assuming PIR
%     dy = H\rhs;
%     return;
%   end;
% end;
%
% if isa(solver, 'function_handle')
%     dy = feval(solver, rhs, H, maxIterCG, tolCG);
%     return
% end

%
%     %     case {'Joint-CG-hyperElastic'}
%     %         Afctn = @(x) H.M(x) + H.d2S.d2S(x,H.omega,H.m,H.d2S.yc);
%     %         [dy,flag,relres,iter] = pcg(Afctn,rhs,tolCG,maxIterCG);
%     %     case {'CG-hyperElastic'}
%     %         M =@(x) H.d2D.P((H.d2D.dr'*H.d2D.d2psi*H.d2D.dr)*H.d2D.P(x));
%     %         Afctn = @(x) M(x) + H.d2S.d2S(x,H.omega,H.m,H.d2S.yc);
%     %         [dy,flag,relres,iter] = pcg(Afctn,rhs,tolCG,maxIterCG);
%     %     otherwise
%     %         if isnumeric(H)
%     %             % if H is a matrix, solve the linear system using MATLAB's backslash
%     %             dy = H\rhs;
%     %         else
%     %             error(solver)
%     %         end
%   otherwise,
%     keyboard
%     error(1)
%     %
%   case {'PCG-hyperElastic'}
%     M =@(x) H.d2D.P((H.d2D.dr'*H.d2D.d2psi*H.d2D.dr)*H.d2D.P(x));
%     Afctn = @(x) M(x) + H.d2S.d2S(x,H.omega,H.m,H.d2S.yc);
%     Ddiag = diag(H.d2D.dr'*H.d2D.d2psi*H.d2D.dr);
%     D = H.d2D.P(full(Ddiag))  +  H.d2S.diag(H.d2S.yc);
%     PC = @(x) D.\x; % Jacobi preconditioner
%     [dy,flag,relres,iter] = pcg(Afctn,rhs,tolCG,maxIterCG,PC);
%
%     %   case {'PCG-elastic','PCG-curvature'}
%     %     % operator
%     %     A         = @(x) ...
%     %       H.d2D.P((H.d2D.dr'*H.d2D.d2psi*H.d2D.dr)*H.d2D.P(x)) ...
%     %       + H.d2S.d2S(x,H.omega,H.m);
%     %     % preconditioner
%     %     Ddiag = diag(H.d2D.dr'*H.d2D.d2psi*H.d2D.dr);
%     %     D = H.d2D.P(full(Ddiag))  +  H.d2S.diag(H.omega,H.m);
%     %     PC = @(x) D.\x; % Jacobi preconditioner
%     %     [dy,flag,relres,iter] = pcg(A,rhs,tolCG,maxIterCG,PC);
%
%
%     %     case {'Joint-CG-hyperElastic'}
%     %         Afctn = @(x) H.M(x) + H.d2S.d2S(x,H.omega,H.m,H.d2S.yc);
%     %         [dy,flag,relres,iter] = pcg(Afctn,rhs,tolCG,maxIterCG);
%     %     case {'CG-hyperElastic'}
%     %         M =@(x) H.d2D.P((H.d2D.dr'*H.d2D.d2psi*H.d2D.dr)*H.d2D.P(x));
%     %         Afctn = @(x) M(x) + H.d2S.d2S(x,H.omega,H.m,H.d2S.yc);
%     %         [dy,flag,relres,iter] = pcg(Afctn,rhs,tolCG,maxIterCG);
%     %     otherwise
%     %         if isnumeric(H)
%     %             % if H is a matrix, solve the linear system using MATLAB's backslash
%     %             dy = H\rhs;
%     %         else
%     %             error(solver)
%     %         end

