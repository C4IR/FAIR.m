%==============================================================================
%
% function [yc,His] = lBFGS(fctn,yc,varargin)
%
% limited BFGS optimizer
%
% Input:
% ------
%   fctn           function handle
%   yc          starting guess 
%   varargin    optional parameter, see below
%
% Output:
% -------
%   yc          numerical optimizer (current iterate)
%   His         iteration history
%==============================================================================

function [yc,His] = lBFGS(fctn,yc,varargin)

if nargin ==0, % help and minimal example
  help(mfilename);  
  E9_MRIhead_MLIRlBFGS_NGF_mbElas;  
  yc = 'endOfMinimalExample'; 
  His = [];
  return;
end;

% parameter initialization -----------------------------------------------
maxIter      = 10;              % maximum number of iterations
tolJ         = 1e-3;            % for stopping, objective function
tolY         = 1e-2;            %   - " -     , current value
tolG         = 1e-2;            %   - " -     , norm of gradient
solver       = regularizer;
maxLBFGS     = 5;              % maximum number of BFGS vectors
LSmaxIter    = 10;              % maximum number of line search iterations
LSreduction  = 1e-4;            % minimal reduction in line search
vecNorm      = @norm;           % norm to be used for dJ and dy    
yStop        = [];              % used for stopping in multi-level framework
Jstop        = [];              % 
Plots        = @(iter,para) []; % for plots;

for k=1:2:length(varargin),     % overwrites default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;
objFctn = fctn([]);

if isempty(yStop), yStop  = yc; end; % yStop used for stopping only
% -- end parameter set-up   ----------------------------------------------

% some output
FAIRmessage([mfilename '(JM 2009/01/31)']);
fprintf('[ maxIter=%d / maxLBFGS=%d / tolJ=%s / tolY=%s / tolG=%s / length(yc)=%d ]\n',...
  maxIter,maxLBFGS,num2str(tolJ),num2str(tolY),num2str(tolG),length(yc));

% -- initialize  ---------------------------------------------------------                                        
STOP = zeros(5,1);

zBFGS = [];     % memory for BFGS gradient directions
sBFGS = [];     % memory for BFGS directions
nBFGS  = 0;     % counter for number of limited BFGS directions

if isempty(Jstop),
  % evaluate objective function for stopping values and plots
  [Jstop,para] = fctn(yStop); Jstop = abs(Jstop) + (Jstop == 0);
  Plots('stop',para);
end;

% evaluate objective function for starting values and plots
[Jc,para,dJ,H0] = fctn(yc); 
Plots('start',para);
iter = 0; yOld = 0*yc; Jold = Jc; y0 = yc;

hisStr    = {'iter','J','Jold-J','|\nabla J|','|dy|','LS'};
his        = zeros(maxIter+2,6);
his(1,1:3) = [-1,Jstop,Jstop-Jc];
his(2,:)   = [0,Jc,Jstop-Jc,vecNorm(dJ),vecNorm(yc-yStop),0];

% some output
fprintf('%4s %-12s %-12s %-12s %-12s %4s\n%s\n',...
  hisStr{:},char(ones(1,64)*'-'));
dispHis = @(var) ...
  fprintf('%4d %-12.4e %-12.3e %-12.3e %-12.3e %4d\n',var);
dispHis(his(1,:));
% -- end initialization   ------------------------------------------------


% ==============================================================================
% MAIN LOOP
% ==============================================================================
while 1,
  % check stopping rules
  STOP(1) = abs(Jold-Jc)  <= tolJ*(1+abs(Jstop));
  STOP(2) = (iter>0) && (norm(yc-yOld) <= tolY*(1+norm(y0)));
  STOP(3) = norm(dJ)      <= tolG*(1+abs(Jstop));
  STOP(4) = norm(dJ)      <= 1e3*eps;
  STOP(5) = (iter >= maxIter);
  if all(STOP(1:3)) || any(STOP(4:5)), break;  end;

  iter      = iter + 1;
  
  
  % update at most maxLBFGS BFGS directions
  if iter > 1,
    zz = (dJ - dJold)';
    ss =  yc - yOld;
    if zz'*ss > 0,      
      start = 2-(nBFGS<maxLBFGS);
      zBFGS = [zBFGS(:,start:end),zz];
      sBFGS = [sBFGS(:,start:end),ss];
      nBFGS = min(maxLBFGS,nBFGS+1);      
    else
      warning('   y''*s < 0, skip update');
    end;
  end;
  
  % compute search direction using recursive and limited BFGS
  dy = bfgsrec(solver,nBFGS,sBFGS,zBFGS,H0,-dJ');
  if iter == 1,
    if strcmp(objFctn,'NPIR'),
      m     = para.m; 
      omega = para.omega; 
      dim   = length(para.omega)/2;
      maxdy = norm(sqrt(sum(reshape(center(dy,m),[],dim).^2,2)),'inf');
      fac   = 0.5*max((omega(2:2:end)-omega(1:2:end)));
      if maxdy>fac,
        dy = (fac/maxdy)*dy;
        %H0 = fac*H0;%speye(length(dy),length(dy));
      end;
    else
      H0 = norm(dJ) * speye(length(dy),length(dy));
    end;
  end;
  % perform Armijo line-search
  [t,yt,LSiter] = Armijo(fctn,yc,dy,Jc,dJ,...
        'LSmaxIter',LSmaxIter,'LSreduction',LSreduction);
  if (t == 0), 
    break; 
  end; % break if line-search fails

  % update variables
  yOld = yc; Jold = Jc; dJold = dJ; yc = yt;
  [Jc,para,dJ] = fctn(yc);  % evaluate objective function
  
  % some output
  his(iter+2,:) = [iter,Jc,Jold-Jc,vecNorm(dJ),vecNorm(yc-yOld),LSiter];
  dispHis(his(iter+1,:));
  para.normdY = vecNorm(yc - yOld);
  Plots(iter,para);
  
end;%while; % end of iteration loop
% ==============================================================================

% clean up
His.str = hisStr;
His.his = his(1:iter+2,:);
fprintf('STOPPING:\n');
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e]\n',STOP(1),...
  '(Jold-Jc)',(Jold-Jc),'tolJ*(1+|Jstop|)',tolJ*(1+abs(Jstop)));
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e]\n',STOP(2),...
  '|yc-yOld|',norm(yc-yOld),'tolY*(1+norm(yc)) ',tolY*(1+norm(yc)));
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e]\n',STOP(3),...
  '|dJ|',norm(dJ),'tolG*(1+abs(Jstop)',tolG*(1+abs(Jstop)));
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e]\n',STOP(4),...
  'norm(dJ)',norm(dJ),'eps',1e3*eps);
fprintf('%d[ %-10s=  %-14d >= %-25s=  %-14d]\n',STOP(5),...
  'iter',iter,'maxIter',maxIter);

FAIRmessage([mfilename,' : done !']);

%------------------------------------------------------------------------------

function[d] = bfgsrec(solver,n,S,Z,H,d)

maxIterCG = 10; tolCG = 1e-2;
if isempty(solver) 
  if isstruct(H),
    if isfield(H,'solver'), 
      solver = H.solver;
    elseif isfield(H,'d2S') && isfield(H.d2S,'solver'), 
      solver = H.d2S.solver;
    else
      error('solver has not been defined')
    end;

  end
end

if n == 0, 
  switch solver,
    case {'MG-elastic'}
        d = MGsolver(d,H);
    case {'pcg'}
      L   = tril(H); % Symmetric Gauss Seidel Preconditioning,
      D   = diag(H); % L is lower, D is diagonal, U = L'
      SGS = @(x) L\(D.*(L'\x));
      d   = pcg(H,d,tolCG,maxIterCG,SGS);
    case 'PCG-curvature',
      M     = @(x) H.d2D.M;
      Afctn = @(x) M(x) + H.d2S.d2S(x,H.omega,H.m);
      % preconditioner
      Ddiag = H.d2D.M;
      D = Ddiag  +  H.d2S.diag(H.omega,H.m);
      PC = @(x) D.\x; % Jacobi preconditioner
      [d,flag,relres,iter] = pcg(Afctn,d,tolCG,maxIterCG,PC);

    otherwise,
      if isnumeric(H),
        % if H is a matrix, solve the linear system using MATLAB's backslash
        d = H\d;
      else
        error('nyi - solver %s', solver)
      end;
  end;   
else
  alpha = (S(:,n)'*d)/(Z(:,n)'*S(:,n));
  d     = d - alpha*Z(:,n);
  d     = bfgsrec(solver,n-1,S,Z,H,d);
  d     = d + (alpha - (Z(:,n)'*d)/(Z(:,n)'*S(:,n)))*S(:,n);
end;
%==============================================================================
