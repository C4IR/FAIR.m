%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function [yc,His] = GaussNewton(fctn,yc,varargin)
%
% Gauss-Newton scheme with variable line search for minimizing J = fctn(yc)
% 
% Input:
%   fctn        function handle
%   yc          starting guess 
%   varargin    optional parameter, see below
%
% Output:
%   yc          numerical optimizer (current iterate)
%   his         iteration history
%
% Note: the linear systems are solved using solveLinearSystem.m
%==============================================================================

function [yc,His] = GaussNewton(fctn,yc,varargin)

if nargin==0
    help(mfilename)
    E9_Hands_NPIR
    yc = 'endOfMinimalExample'; 
    return;
end;

% parameter initialization -----------------------------------------------
maxIter      = 10;              % maximum number of iterations
tolJ         = 1e-3;            % for stopping, objective function
tolY         = 1e-2;            %   - " -     , current value
tolG         = 1e-2;            %   - " -     , norm of gradient
lineSearch   = @Armijo;         % linesearch scheme
LSmaxIter    = 10;              % maximum number of line search iterations
LSreduction  = 1e-4;            % minimal reduction in line search
vecNorm      = @norm;           % norm to be used for dJ and dy    
solver       = [];              % linear solver 
yStop        = [];              % used for stopping in multi-level framework
Jstop        = [];              % 
Plots        = @(iter,para) []; % for plots;

for k=1:2:length(varargin),     % overwrites default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

if ~isa(Plots,'function_handle') && (Plots == 0 || strcmp(Plots,'off')),
  Plots        = @(iter,para) []; % for plots;
end;

if isempty(yStop), yStop  = yc; end; % yStop used for stopping only
% -- end parameter setup   ----------------------------------------------

% some output
if isstring(solver)
    solverName = solver;
elseif isa(solver,'function_handle')
    solverName = func2str(solver);
else
    error(1);
end
FAIRmessage(sprintf('JM 2011/01/02 : %s / %s / %s',...
  mfilename,func2str(lineSearch),solverName))
fprintf('[ maxIter=%s / tolJ=%s / tolY=%s / tolG=%s / numel(yc)=%d]\n',...
  num2str(maxIter),num2str(tolJ),num2str(tolY),num2str(tolG),length(yc));

% -- initialize  ---------------------------------------------------------                                        
STOP = zeros(5,1);

if isempty(Jstop),
  % evaluate objective function for stopping values and plots
  [Jstop,para] = fctn(yStop); Jstop = abs(Jstop) + (Jstop == 0);
  Plots('stop',para);
end;

% evaluate objective function for starting values and plots
[Jc,para,dJ,H] = fctn(yc); 
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
dispHis(his(2,:));
% -- end initialization   ------------------------------------------------


%-- start the iteration --------------------------------------------------
while 1, 
  % check stopping rules
  STOP(1) = (iter>0) && abs(Jold-Jc)   <= tolJ*(1+abs(Jstop));
  STOP(2) = (iter>0) && (norm(yc-yOld) <= tolY*(1+norm(y0)));
  STOP(3) = norm(dJ)                   <= tolG*(1+abs(Jstop));
  STOP(4) = norm(dJ)                   <= 1e6*eps;
  STOP(5) = (iter >= maxIter);
  if all(STOP(1:3)) || any(STOP(4:5)), break;  end;

  iter = iter + 1;
  
  % solve the Gauss-Newton System
  [dy,solver] = solveLinearSystem(-dJ',H,solver);
  
  % check descent direction
  % note: descent is not granted if using an iterative solver 
  descent =   dJ * dy; 
  if descent > 0,
    warning('no descent direction, switch to -dy!') 
    dy      = -dy;
  end;

  % perform Armijo line-search
  [t,yt,LSiter] = lineSearch(fctn,yc,dy,Jc,dJ,...
      'para',para,'LSmaxIter',LSmaxIter,'LSreduction',LSreduction);
  if (t == 0),
      break; 
  end; % break if line-search fails
  
  % save old values and update
  yOld = yc; Jold = Jc; yc = yt;
  [Jc,para,dJ,H] = fctn(yc); % evalute objective function
  
  % some output
  his(iter+2,:) = [iter,Jc,Jold-Jc,vecNorm(dJ),vecNorm(yc-yOld),LSiter];
  dispHis(his(iter+2,:));
  para.normdY = vecNorm(yc - yOld);
  Plots(iter,para);
% pause
end;%while; % end of iteration loop
%-------------------------------------------------------------------------
Plots(iter,para);

% clean up
His.str = hisStr;
His.his = his(1:iter+2,:);
fprintf('STOPPING:\n');
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e]\n',STOP(1),...
  '(Jold-Jc)',(Jold-Jc),'tolJ*(1+|Jstop|)',tolJ*(1+abs(Jstop)));
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e]\n',STOP(2),...
  '|yc-yOld|',norm(yc-yOld),'tolY*(1+norm(yc)) ',tolY*(1+norm(yc)));
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e]\n',STOP(3),...
  '|dJ|',norm(dJ),'tolG*(1+abs(Jstop))',tolG*(1+abs(Jstop)));
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e]\n',STOP(4),...
  'norm(dJ)',norm(dJ),'eps',1e3*eps);
fprintf('%d[ %-10s=  %-14d >= %-25s=  %-14d]\n',STOP(5),...
  'iter',iter,'maxIter',maxIter);

FAIRmessage([mfilename,' : done !'],'=');

%==============================================================================
