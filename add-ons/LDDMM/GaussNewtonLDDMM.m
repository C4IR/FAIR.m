%==============================================================================
% This code is part of the Matlab-based toolbox LagLDDDM - A Lagrangian Gauss--
% Newton--Krylov Solver for Mass- and Intensity-Preserving Diffeomorphic Image
% Registration
%
% For details and license info see
% - https://github.com/C4IR/LagLDDMM%
%==============================================================================
%
% function [yc,His] = GaussNewtonLDDMM(fctn,yc,varargin)
%
% Gauss-Newton scheme tailored for LDDMM. The numerics is the same as in
% GaussNewton.m, but the stopping criterion and some outputs are modified.
%
% Input:
%   fctn        function handle
%   vc          initial guess
%   varargin    optional parameter, see below
%
% Output:
%   vc          numerical optimizer (current iterate)
%   his         iteration history
%
% see also GaussNewton.m
%==============================================================================

function [vc,His] = GaussNewtonLDDMM(fctn,vc,varargin)

if nargin ==0, % help and minimal example
  help(mfilename);  E10_2Ddisc2C_hyperElastic;  vc = 'example finished'; return;
end;

% parameter initialization -----------------------------------------------
lineSearch   = @Armijo;         % default line search
maxIter      = 10;              % maximum number of iterations
tolJ         = 1e-3;            % for stopping, objective function
tolY         = 1e-2;            %   - " -     , current value
tolG         = 1e-1;            %   - " -     , norm of gradient
LSMaxIter    = 10;              % maximum number of line search iterations
LSreduction  = 1e-4;            % minimal reduction in line search
vecNorm      = @norm;           % norm to be used for dJ and dy
solver       = [];              % linear solver
yStop        = [];
vStop        = [];              % used for stopping in multi-level framework
Jstop        = [];              %
Plots        = @(iter,para) []; % for plots;
for k=1:2:length(varargin),     % overwrites default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

%fprintf = @(varargin) [];

if ~isa(Plots,'function_handle') && (Plots == 0 || strcmp(Plots,'off')),
  Plots        = @(iter,para) []; % for plots;
end;

if isempty(vStop), vStop  = vc; end; % vStop used for stopping only
% -- end parameter setup   ----------------------------------------------

% some output
% FAIRmessage = @(str) fprintf('%% %s  [ %s ]  % s\n',...
%   char(ones(1,10)*'-'),str,char(ones(1,60-length(str))*'-'));
FAIRmessage([mfilename '(LR 2017/06/01)']);
fprintf('[ maxIter=%s / tolJ=%s / tolY=%s / tolG=%s / length(yc)=%d ]\n',...
  num2str(maxIter),num2str(tolJ),num2str(tolY),num2str(tolG),length(vc));

% -- initialize  ---------------------------------------------------------
STOP = zeros(5,1);

if isempty(Jstop) || isempty(yStop),
  % evaluate objective function for stopping values and plots
  [Jstop,para] = fctn(yStop); Jstop = abs(Jstop) + (Jstop == 0); yStop = para.yc;
  Plots('stop',para);
end;

% evaluate objective function for starting values and plots
[Jc,para,dJ,H] = fctn(vc);
Jac = geometry(para.yc,para.m-1,'Jac','matrixFree',true,'omega',para.omega);
Plots('start',para);
iter = 0; yOld = 0*vc; Jold = Jc; y0 = para.yc;

hisStr    = {'iter','J','Jold-J','|\nabla J|','|dy|','LS','min(Jac)','max(Jac)','nnz(H)','iterCG','relres'};
his        = zeros(maxIter+2,11);
his(1,1:3) = [-1,Jstop,Jstop-Jc];

nnzH = -1;
if isnumeric(H); nnzH = nnz(H); end;
his(2,:)   = [0,Jc,Jstop-Jc,vecNorm(dJ),vecNorm(y0-yStop),0,min(Jac),max(Jac),nnzH,0,0.0];

% some output
fprintf('%4s %-12s %-12s %-12s %-12s %4s  %-8s  %-8s  %-5s %-12s%-5s\n%s\n',...
  hisStr{:},char(ones(1,111)*'-'));
dispHis = @(var) ...
  fprintf('%4d %-12.4e %-12.3e %-12.3e %-12.3e %-4d %-1.4f\t%-1.4f\t  %-12d %-3d   %-12.3e\n',var);
dispHis(his(1,:));
dispHis(his(2,:));
% -- end initialization   ------------------------------------------------


%-- start the iteration --------------------------------------------------
while 1,
  % check stopping rules
  STOP(1) = (iter>0) && abs(Jold-Jc)   <= tolJ*(1+abs(Jstop));
  STOP(2) = (iter>0) && (norm(para.yc-yOld) <= tolY*(1+norm(y0)));
  STOP(3) = norm(dJ)                   <= tolG*(1+abs(Jstop));
  STOP(4) = norm(dJ)                   <= 1e6*eps;
  STOP(5) = (iter >= maxIter);
  if all(STOP(1:3)) || any(STOP(4:5)), break;  end;

  iter = iter + 1;
  % solve the Gauss-Newton System
  [dy,it,rr] = solveGN(-dJ',H,solver);

  % check descent direction
  % note: descent is not granted if using an iterative solver
  descent =   dJ * dy;
  if descent > 0,
    warning('no descent direction, switch to -dy!') %#ok<WNTAG>
    dy      = -dy;
  end;

  % perform Armijo line-search
  [t,yt,LSiter] = lineSearch(fctn,vc,dy,Jc,dJ,'para',para,...
        'LSMaxIter',LSMaxIter,'LSreduction',LSreduction);


  if LSiter>2
%        keyboard;
  end
  if (t == 0),
      break;
  end; % break if line-search fails

  % save old values and update
  yOld = para.yc; Jold = Jc; vc = yt;
  [Jc,para,dJ,H] = fctn(vc); % evalute objective function
  Jac = geometry(para.yc,para.m-1,'Jac','matrixFree',true,'omega',para.omega);

%   clf; plot(para.T,para.d,'r',para.T,para.g,'g')
%   pause

  % some output
  if isnumeric(H); nnzH = nnz(H); end;
  his(iter+2,:) = [iter,Jc,Jold-Jc,vecNorm(dJ),vecNorm(para.yc-yOld),LSiter,min(Jac),max(Jac),nnzH,it,rr];
  dispHis(his(iter+2,:));
  para.normdY = vecNorm(para.yc - yOld);
  Plots(iter,para);
% pause
end;%while; % end of iteration loop
%-------------------------------------------------------------------------
Plots(iter,para);

% clean up
His.str = hisStr;
His.his = his(1:iter+2,:);
His.para = para;
fprintf('STOPPING:\n');
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e]\n',STOP(1),...
  '(Jold-Jc)',(Jold-Jc),'tolJ*(1+|Jstop|)',tolJ*(1+abs(Jstop)));
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e]\n',STOP(2),...
  '|yc-yOld|',norm(para.yc-yOld),'tolY*(1+norm(yc)) ',tolY*(1+norm(y0)));
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e]\n',STOP(3),...
  '|dJ|',norm(dJ),'tolG*(1+abs(Jstop))',tolG*(1+abs(Jstop)));
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e]\n',STOP(4),...
  'norm(dJ)',norm(dJ),'eps',1e3*eps);
fprintf('%d[ %-10s=  %-14d >= %-25s=  %-14d]\n',STOP(5),...
  'iter',iter,'maxIter',maxIter);

FAIRmessage([mfilename,' : done !']);

%==============================================================================

function [dy,iter,relres] = solveGN(rhs,H,solver)
maxIterCG = 50; tolCG = 1e-1;
iter = -1; relres = 0.0;
if isempty(solver)
  if isstruct(H),
    if isfield(H,'solver'),
      solver = H.solver;
    elseif isfield(H,'d2S') && isfield(H.d2S,'solver'),
      solver = H.d2S.solver;
    else
      error('solver has not been defined')
    end;
  else % no regularizer initialized, assuming PIR
    dy = H\rhs;
    return;
  end;
end;

if isa(solver, 'function_handle')
    [dy,iter,relres] = feval(solver, rhs, H, maxIterCG, tolCG);
    return
end

switch solver
    % matrix based
    % ------------
    case 'pcg'
        L   = tril(H); % Symmetric Gauss Seidel Preconditioning,
        D   = diag(H); % L is lower, D is diagonal, U = L'
        SGS = @(x) L\(D.*(L'\x));
        [dy,flag,relres,iter,resvec] = pcg(H,rhs,tolCG,maxIterCG,SGS);
    case 'ichol'
        L1   = cholinc(sparse(H),'0'); % Symmetric Gauss Seidel Preconditioning,
        [dy,flag,relres,iter,resvec] = pcg(H,rhs,tolCG,maxIterCG,L1,L1');
    case 'jacobi-pcg'
        D   = diag(H); % D is diagonal
        PC = @(x) D.\x;
        [dy,flag,relres,iter,resvec] = pcg(H,rhs,tolCG,maxIterCG,PC);
    case 'cg'
        [dy,flag,relres,iter,resvec] = pcg(H,rhs,tolCG,maxIterCG);
    case 'PIR-cg'
        if isstruct(H)
            [dy,flag,relres,iter,resvec] = pcg(H.operator,rhs,tolCG,maxIterCG);
        else
            [dy,flag,relres,iter,resvec] = pcg(H,rhs,tolCG,maxIterCG);
        end
    % matrix free
    % ------------
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
    % whereas in conjugate gradient (CG, available for elastic, curvature and hyperelastic) methods we can
    % model also the off-diagonals d2D(x) = (P'*dr'*d2psi*dr*P) * x
    %
    % The default solver for a chosen regularizer is parameterized by d2S.solver.
    % However, one can overload this setting using
    %    >> MLIR(..., 'solverNPIR','myFavoriteSolver', ... );
    case {'MG-elastic'}
        [dy,relres,iter] = MGsolver(rhs,H);
    case {'CG-elastic','CG-curvature'}
        M         = @(x) H.d2D.P((H.d2D.dr'*H.d2D.d2psi*H.d2D.dr)*H.d2D.P(x));
        Afctn     = @(x) M(x) + H.d2S.d2S(x,H.omega,H.m);
        [dy,flag,relres,iter,resvec] = pcg(Afctn,rhs,tolCG,maxIterCG);
    case {'PCG-elastic','PCG-curvature'}
        M         = @(x) H.d2D.P(H.d2D.dr'*(H.d2D.d2psi*(H.d2D.dr*H.d2D.P(x))));
        Afctn     = @(x) M(x) + H.d2S.d2S(x,H.omega,H.m);
        % preconditioner
%         Ddiag = diag(H.d2D.dr'*H.d2D.d2psi*H.d2D.dr);
        Ddiag = sum(H.d2D.d2psi*H.d2D.dr.^2,1);
        Ddiag = Ddiag(:);
        D = H.d2D.P(full(Ddiag))  +  H.d2S.diag(H.omega,H.m);
        PC = @(x) D.\x; % Jacobi preconditioner
        [dy,flag,relres,iter,resvec] = pcg(Afctn,rhs,tolCG,maxIterCG,PC);
    case {'Joint-CG-hyperElastic'}
        Afctn = @(x) H.M(x) + H.d2S.d2S(x,H.omega,H.m,H.d2S.yc);
        [dy,flag,relres,iter,resvec] = pcg(Afctn,rhs ,tolCG,maxIterCG);
    case {'CG-hyperElastic'}
        M =@(x) H.d2D.P((H.d2D.dr'*H.d2D.d2psi*H.d2D.dr)*H.d2D.P(x));
        Afctn = @(x) M(x) + H.d2S.d2S(x,H.omega,H.m,H.d2S.yc);
        [dy,flag,relres,iter,resvec] = pcg(Afctn,rhs,tolCG,maxIterCG);
    case {'PCG-hyperElastic'}
        M =@(x) H.d2D.P((H.d2D.dr'*H.d2D.d2psi*H.d2D.dr)*H.d2D.P(x));
        Afctn = @(x) M(x) + H.d2S.d2S(x,H.omega,H.m,H.d2S.yc);
        Ddiag = diag(H.d2D.dr'*H.d2D.d2psi*H.d2D.dr);
        D = H.d2D.P(full(Ddiag))  +  H.d2S.diag(H.d2S.yc);
        PC = @(x) D.\x; % Jacobi preconditioner
        [dy,flag,relres,iter,resvec] = pcg(Afctn,rhs,tolCG,maxIterCG,PC);
    otherwise
        if isnumeric(H)
            % if H is a matrix, solve the linear system using MATLAB's backslash
            dy = H\rhs;
        else
            error(solver)
        end
end

if exist('flag','var')
    switch flag
        case 1
            fprintf('pcg iterated %d times without converging to tolerance %e. Returned iterate (number %d) has residual %e\n',...
                numel(resvec)-1,tolCG,iter,relres);
            iter = numel(resvec)-1;
        case 2
            fprintf('Preconditioner M was ill-conditioned.\n');
        case 3
            fprintf('pcg stagnated. (Two consecutive iterates were the same.)\n');
        case 4
            fprintf('One of the scalar quantities calculated during pcg became too small or too large to continue computing.\n');
        otherwise
%             fprintf('pcg success! %d iterations / relres= %1.2e / tolCG= %1.2e\n',iter,relres,tolCG);
    end
end
