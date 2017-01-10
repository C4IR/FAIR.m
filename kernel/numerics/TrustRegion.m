%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function [Yc,His] = TrustRegion(fctn,Yc,varargin)
%
% Trust-Region scheme for minimizing J = fctn(Yc)
% 
% Input:
%   fctn         function handle
%  Yc      starting guess 
%   varargin    optional parameter, see below
%
% Output:
%   Yopt      numerical optimizer
%   His         iteration history
%==============================================================================

function [Yc,His] = TrustRegion(fctn,Yc,varargin)

if nargin ==0, 
    help(mfilename); 
    E9_Hands_NPIRmf_TR_nopre; 
    Yc = 'endOfMinimalExample';
    return; 
end;

maxIter   = 10;             % maximum number of iterations
tolJ      = 1e-2;           % for stopping, objective function
tolY      = 5e-3;           %   - " -     , current value
tolG      = 1e-0;           %   - " -     , norm of gradient
Ystop     = [];             % used for stopping in multi-level framework
vecNorm   = @norm;          % norm for transformations
preconditioner = @(x,para) x;    % Linear solver: precondtioner 
LSmaxIter = 10;             % Linear solver: maximum number of iterations
LStol     = 0.1;            % Linear solver:  tolerance
TRmaxIter = 10;             % TR: maximum number of iterations
TRmu      = .25;            % TR: success if ratio > TRmu
TRshrink  = .25;            % TR: shrink radius if ratio <  TRenlarge
TRenlarge = .75;            % TR: enlarge radius if ratio >= TRenlarge
enlarge   = 2;              % TR: enlarge trust region by enlarge
shrink    = 0.5;            % TR: shrink  trust region by shrink
Plots   = @(iter,para) [];  % for plots;

for k=1:2:length(varargin), % overwrites default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

if not(regularizer('get','matrixFree')),
    error('%s - matrixBased regularizers nyi.',mfilename);
end
if strcmp(regularizer,'mfHyperElastic'),
    error('%s - hyperElastic regularizer nyi.',mfilename);
end


if isempty(Ystop),  Ystop  = Yc;    end; % Ystop: used for stopping only
% -- end parameter set-up   ----------------------------------------------

% some output
FAIRmessage = @(str) fprintf('%% %s  [ %s ]  % s\n',...
  char(ones(1,10)*'-'),str,char(ones(1,60-length(str))*'-'));
FAIRmessage([mfilename '(JM 2008/08/08)']);
fprintf('/maxIter=%s / tolJ=%s / tolY=%s / tolG=%s / length(Yc)=%d/\n',...
  num2str(maxIter),num2str(tolJ),num2str(tolY),num2str(tolG),length(Yc));

% -- initialize  ----------------------------------------------------------
STOP = zeros(5,1);
TRradius = max(norm(Yc),1);             % Initialize TR size

% evaluate objective function for stopping values and plots
[Jstop,para] = fctn(Ystop); Jstop = abs(Jstop) + (Jstop == 0);
Plots('stop',para);

% evaluate objective function for starting values and plots
[Jc,para,dJ,H] = fctn(Yc);
% set up Hessian of data term (I assume that regularizer is always mf!!!)
H.d2D.M = @(x) H.d2D.P(H.d2D.dr' * H.d2D.d2psi * H.d2D.dr * H.d2D.P(x));

iter = 0; Yold = 0*Yc; Jold = Jc; Y0 = Yc;
Plots('start',para);

His.str = {'iter','J','Jold-J','|\nabla J|','|dY|','TR'};
his    = zeros(maxIter+2,6);
his(1,1:3) = [-1,Jstop,Jstop-Jc];
his(2,:)   = [0,Jc,Jstop-Jc,vecNorm(dJ),vecNorm(Yc-Ystop),0];

% some output
fprintf('%4s %-12s %-12s %-12s %-12s %4s\n%s\n',...
  His.str{:},char(ones(1,64)*'-'));
dispHis = @(var) ...
  fprintf('%4d %-12.4e %-12.3e %-12.3e %-12.3e %4d\n',var);
dispHis(his(1,:));
% -- end initializarion   -------------------------------------------------

%-- start the iteration ---------------------------------------------------
while 1,
  % check stopping rules
  STOP(1) = abs(Jold-Jc)  <= tolJ*(1+abs(Jstop));
  STOP(2) = norm(Yc-Yold) <= tolY*(1+norm(Y0));
  STOP(3) = norm(dJ)      <= tolG*(1+abs(Jstop));
  STOP(4) = norm(dJ)      <= 1e3*eps;
  STOP(5) = (iter > maxIter);
  if all(STOP(1:3)) || any(STOP(4:5)), break;  end;

  iter = iter + 1;
  
  % solve Gauss-Newton system with CG
  [dY,Xcg] = TRcg(H,-dJ',TRradius,LStol,LSmaxIter,preconditioner);
  H.d2D.M = @(x) H.d2D.P(H.d2D.dr' * H.d2D.d2psi * H.d2D.dr * H.d2D.P(x));

%  relNorm = norm(mfAy(dY,H)+dJ')/norm(dJ)
  Yt    = Yc+dY;
  Jt    = fctn(Yt);
  ared  = Jt-Jc;
  pred  = dJ*dY + .5*dY'*mfAy(dY,H);
  ratio = ared/pred;

  if (ratio > TRmu), % step is working
    if (ratio > TRenlarge) && (norm(dY) > (TRradius-1e-8)),
      TRradius = enlarge*TRradius;
    elseif (ratio < TRshrink)
      TRradius = shrink*TRradius;
    end
  else             % problem: step need to be reduced
    TRiter = 1;
    while (ratio <= TRmu) && (TRiter < TRmaxIter)
      TRradius = shrink*min(TRradius,norm(dY));
      dY = TRadjust(TRradius,Xcg);
      Yt = Yc + dY;
      Jt = fctn(Yt);
      ared = Jt - Jc;
      % compute new pred only if ared < 0
      if ared < 0,
        pred = dJ*dY + .5*dY'*mfAy(dY,H);
      end
      ratio  = ared/pred;
      TRiter = TRiter+1;
      if TRiter > TRmaxIter, % big problem: better have a beer
        warning('Trust Region approach failed!')
        return
      end
    end
  end;

  Yold = Yc; Jold = Jc; Yc = Yt;    % save old values and update
  [Jc,para,dJ,H] = fctn(Yc);        % evalute objective function


  % some output
  his(iter+2,:) = [iter,Jc,Jold-Jc,vecNorm(dJ),vecNorm(Yc-Yold),TRradius];
  dispHis(his(iter+1,:));
  para.normdY = vecNorm(Yc - Yold);
  Plots(iter,para);

end;%while; % end of iteration loop
%------------------------------------------------------------------------------

% clean up
His.his = his(1:iter+1,:);

fprintf('STOPPING:\n');
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e]\n',STOP(1),...
  '(Jold-Jc)',(Jold-Jc),'tolJ*(1+|Jstop|)',tolJ*(1+abs(Jstop)));
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e]\n',STOP(2),...
  '|Yc-Yold|',norm(Yc-Yold),'tolY*(1+norm(Yc)) ',tolY*(1+norm(Yc)));
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e]\n',STOP(3),...
  '|dJ|',norm(dJ),'tolG*(1+abs(Jstop)',tolG*(1+abs(Jstop)));
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e]\n',STOP(4),...
  'norm(dJ)',norm(dJ),'eps',1e3*eps);
fprintf('%d[ %-10s=  %-14d >= %-25s=  %-14d]\n',STOP(5),...
  'iter',iter,'maxIter',maxIter);

FAIRmessage([mfilename,' : done !']);

%------------------------------------------------------------------------------

function [s,k] = TRadjust(radius, X)
% find the intersetion s of the path pk = sum(X(:,k)), k=1,...,
% with the trustregion boundary

s = X(:,1);              % start with stepest descent
if (norm(s) > radius) || (size(X,2) == 1),
  s = (radius/norm(s))*s;        % Cauchy Point
  k = 1;
  return;
end;

for k=2:size(X,2),
  if (norm(X(:,k)) <= radius),
    % step can be enlarged
    s = X(:,k);
  else
    % compute intersection |s+alpha*d|^2=radius^2
    d = X(:,k)-s;
    a = d'*d;
    b = 2*(s'*d);
    c = s'*s - radius*radius;
    alpha = (-b + sqrt(b*b - 4*a*c))/(2*a);
    s = s +alpha*d;
    return;
  end;
end;

%------------------------------------------------------------------------------

function [x,X]  = TRcg(H,b,radius,tol,maxIter,precond)
% solves x=H\(-g), ||x||<radius
% Input:    
%   H       matrix
%   g       current gradient
%   radius  TR radius
%   pfun    function handle to preconditioner px = pfun(x), default: px=x
%   prec    precondtioner parameters
%
% Output:   
%   x           trial step
%   directions  2D-array of search directions for TR radius reduction

% initialization
n = length(b);
x = zeros(n,1); r = b;

normr0 = norm(b);
normr  = normr0;
iter   = 0;
X      = zeros(n,maxIter);

while((normr > tol*normr0) && (iter <= maxIter) && norm(x) <= radius*(1-1e-8))
  iter = iter + 1;
  if (iter>1), rhoold = rho; end;
  % Preconditioning
  H.d2D.M = 0;% diagonal of Hessian
  z = precond(r,H); % either identity or multigrid
  
  rho = z'*r;
  if (iter==1)
    p = z;
  else
    beta = rho/rhoold;
    p = z + beta*p;
  end
  H.d2D.M = @(x) H.d2D.P(H.d2D.dr' * H.d2D.d2psi * H.d2D.dr * H.d2D.P(x));
  w = mfAy(p,H);
  pTw = p'*w;

  % If alpha <=0 head to the TR boundary and return

  stop = 0;
  if (pTw <= 0)
    a = p'*p; b = 2*(x'*p); c = x'*x - radius*radius;
    alpha = (-b + sqrt(b*b - 4*a*c))/(2*a);
    warning(' negative curvature ')
    stop = 1;
  else
    alpha = rho/pTw;
    if norm(x+alpha*p) > radius
      a = p'*p; b = 2*(x'*p); c = x'*x - radius*radius;
      alpha = (-b + sqrt(b*b - 4*a*c))/(2*a);
      stop = 1;
    end
  end

  x = x+alpha*p;
  X(:,iter) = x;
  r = r - alpha*w;
  normr = norm(r);
  if stop, break, end;
end
X = X(:,1:iter);

%==============================================================================
