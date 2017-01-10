%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function [Yc,His] = SteepestDescent(fctn,Yc,varargin);
%
% Steepest Descent scheme with line search for minimizing J = fctn(Yc)
% 
% Input:
%   fctn         function handle
%   Yc           starting guess 
%   varargin     optional parameter, see below
%
% Output:
%   Yopt         numerical optimizer
%   His          iteration history
%==============================================================================

function [Yopt,His] = SteepestDescent(fctn,Yc,varargin)

if nargin ==0, 
  help(mfilename); 
  % run minimal Example
  E6_Hands_PIR_SD; 
  Yopt = 'endOfMinimalExample';
  return; 
end;

% -- start parameter setup ----------------------------------------------
Ystop       = [];       % default reference parameter for stopping
stepLength  = 1;        % scaling of gradient direction
vecNorm      = @norm;           % norm to be used for dJ and dy    
tolJ         = 1e-3;            % for stopping, objective function
tolY         = 1e-2;            %   - " -     , current value
tolG         = 1e-2;            %   - " -     , norm of gradient
lineSearch   = @Armijo;
LSmaxIter    = 10;              % maximum number of line search iterations
LSreduction  = 1e-4;            % minimal reduction in line search
Plots      = @(iter,para) []; % for plots;

for k=1:2:length(varargin), % overwrites default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

if isempty(Ystop), Ystop = Yc;          end;    % check Ystop
his    = zeros(maxIter+2,6);                    % summarizes iteration history

FAIRmessage(sprintf('JM 2008/08/06 : %s / %s ',mfilename,func2str(lineSearch)))
% report objective function status and prepare plots
fctn([]);
fprintf('  %20s : %s\n','maxIter',num2str(maxIter));
fprintf('\n\n')
% -- end parameter setup   ----------------------------------------------


% -- start initial phase -------------------------------------------------
STOP = zeros(5,1);
[Jstop,para]  = fctn(Ystop);            % compute reference for stopping
Plots('stop',para);

iter = 0; Jold = 0; Yold = 0;           % set initial values
[Jc,para,dJ]  = fctn(Yc);               % compute current values
Plots('start',para);

his(1,:)  = [-1,Jstop,0,0,0,0]; 
his(2,:)  = [0,Jc,abs(Jc-Jstop),vecNorm(dJ),vecNorm(Yc-Yold),0]; 

% headlines for the output table
fprintf('%-4s %-12s %-12s %-12s %-12s %-4s\n',...
  'iter','J','Jold-Jc','|dJ|','|dy|','LS');
fprintf('%-4d %-12.2e %-12.2e %-12.2e %-12.2e %-4d\n',his(1,:));
fprintf('%-4d %-12.2e %-12.2e %-12.2e %-12.2e %-4d\n',his(2,:));
% -- end initial phase   -------------------------------------------------

t = 1; % initial step
if stepLength == 1, stepLength = 0.5/norm(dJ); end;
% -- start iteration phase -----------------------------------------------
while 1,
  tOld = t;
  
  % check stopping rules
  STOP(1) = abs(Jold-Jc)  <= tolJ*(1+abs(Jstop));
  STOP(2) = norm(Yc-Yold) <= tolY*(1+norm(Yc));
  STOP(3) = norm(dJ)      <= tolG*(1+abs(Jstop));
  STOP(4) = norm(dJ)      <= 1e3*eps;
  STOP(5) = (iter >= maxIter);
  if all(STOP(1:3)) | any(STOP(4:5)), break;  end;
  
  iter    = iter + 1;
  
  % scaled steepest descent direction
  dY  = -stepLength*dJ';
  
  % perform line-search

  [t,Yt,LSiter] = lineSearch(fctn,Yc,dY,Jc,dJ,....
    'LSmaxIter',LSmaxIter,'LSreduction',LSreduction);
  if (t == 0), break; end; % break if line-search fails
  
  if (t == 1),
    % everything is fine, be more optimistic
    stepLength = 2*stepLength;  
  else
    % adjust stepLength for the next step
    stepLength = stepLength*t;  
  end;
  
  % update   
  Yold = Yc; Jold = Jc;  Yc = Yt;
  [Jc,para,dJ] = fctn(Yc);
  para.normdY = vecNorm(Yc-Yold);
  Plots(iter,para);

  % store intermediate results, some output
  his(iter+2,:)  = [iter,Jc,abs(Jc-Jold),norm(dJ),para.normdY,LSiter]; 
  fprintf('%-4d %-12.2e %-12.2e %-12.2e %-12.2e %-4d\n',his(iter+2,:));
  
end;%while
% -- end iteration phase -------------------------------------------------

Yopt = Yc;
His.his = his(1:iter+1,:);
His.str = {'iter','J','Jold-Jc','|\nabla J|','|dY|','LS'};

fprintf('\n');
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e ]\n',STOP(1),...
  '(Jold-Jc)',(Jold-Jc),'tolJ*(1+|Jstop|)',tolJ*(1+abs(Jstop)));
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e ]\n',STOP(2),...
  '|Yc-Yold|',norm(Yc-Yold),'tolY*(1+norm(Yc)) ',tolY*(1+norm(Yc)));
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e ]\n',STOP(3),...
  '|dJ|',norm(dJ),'tolG*(1+abs(Jstop)',tolG*(1+abs(Jstop)));
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e ]\n',STOP(4),...
  'norm(dJ)',norm(dJ),'eps',1e3*eps);
fprintf('%d[ %-10s=%-16d  >= %-25s=%-16d ]\n',STOP(5),...
  'iter',iter,'maxIter',maxIter);

FAIRmessage([mfilename,' : done !']);
%==============================================================================
