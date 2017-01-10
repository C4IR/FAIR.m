%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% function [Yopt,His] = Nesterov(fctn,Yc,varargin)
%
% Nesterov-scheme for minimizing J = fctn(yc)
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
%==============================================================================

function [Yopt,His] = Nesterov(fctn,Yc,varargin)

if nargin ==0, 
  help(mfilename); 
  % run minimal Example
  E6_Hands_PIR_Nesterov; 
  Yopt = 'endOfMinimalExample';
  return; 
end;

% % -- start parameter setup ----------------------------------------------
yStop       = [];       % default reference parameter for stopping
maxIter     = 20;       % max number of iterations
tolJ        = 1e-16;     % stopping: relative variation of objective function
tolG        = 5e-8;     % stopping: relative variation of the gradient
tolY        = 5e-8;     % stopping: relative variation of the parameters
vecNorm     = @norm;
LSmaxIter   = 10;       % LineSearch: max number of line search steps
lipschitzC  = 0;
% LSreduction = 1e-6;     % LineSearch: guaranteed reduction by line search
stepLength  = 1;        % scaling of gradient direction
Plots       = @(iter,para) []; % for plots;
Plots       = @FAIRplots;
for k=1:2:length(varargin), % overwrites default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

if isempty(yStop), yStop = Yc;          end;    % check Ystop
his    = zeros(maxIter+2,6);                    % summarizes iteration history

FAIRmessage(sprintf('JM 2008/08/06 : %s ',mfilename))

% report objective function status and prepare plots
fctn([]);
fprintf('  %20s : %s\n','maxIter',num2str(maxIter));
fprintf('\n\n')
% -- end parameter setup   ----------------------------------------------


% -- start initial phase -------------------------------------------------
STOP = zeros(5,1);
[Jstop,para]  = fctn(yStop);            % compute reference for stopping
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
% -- end initial phase   -------------------------------------------------

Yold  = Yc;
Ynold = Yc;
if lipschitzC == 0,
  lipschitzC = norm(dJ)/2
end;

% -- start iteration phase -----------------------------------------------
while 1,
  iter    = iter + 1;
  
%   % check stopping rules
  STOP(1) = abs(Jold-Jc)  <= tolJ*(1+abs(Jstop));
  STOP(2) = norm(Yc-Yold) <= tolY*(1+norm(Yc));
  STOP(3) = norm(dJ)      <= tolG*(1+abs(Jstop));
  STOP(4) = norm(dJ)      <= 1e3*eps;
  STOP(5) = (iter > maxIter);
  if all(STOP(1:3)) | any(STOP(4:5)), break;  end;
  
%   % Nesterov 1 
%   Yn  = Yc + (iter-1)/(iter+2)*(Yc -Yold);

  % Kim & Fessler 2 
  Yn  = Yc + (iter-1)/(iter+2)*(Yc -Yold) + (iter+1)/(iter+2) * (Yc-Ynold);
  [Jn,para,dJ]  = fctn(Yn);

  OK = 0;
  metric = @(Yt) norm(Yt-Yn)^2; 
  for k=1:LSmaxIter,
      Yt = Yn - 1/lipschitzC * dJ';
      Jt = fctn(Yt);
      
      if Jt <= Jn + dJ*(Yt - Yn) + lipschitzC/2 * metric(Yt);
          OK = 1;
          lipschitzC = lipschitzC/2;
          break
      else
          lipschitzC = 2*lipschitzC;
      end
  end
  
  if ~OK,
      disp('LS error');
      keyboard
  end;
  
  % update   
  Ynold = Yn; 
  Yold  = Yc;   Jold = Jc;  
  Yc = Yt; Jc = Jt;
  para.normdY = vecNorm(Yc-Yold);
  Plots(iter,para);

  % store intermediate results, some output
  his(iter+2,:)  = [iter,Jc,abs(Jc-Jold),norm(dJ),para.normdY,lipschitzC]; 
  fprintf('%-4d %-12.2e %-12.2e %-12.2e %-12.2e %-4d\n',his(iter+1,:));
  
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
fprintf('%d[ %-10s=%-16d  > %-25s=%-16d ]\n',STOP(5),...
  'iter',iter,'maxIter',maxIter);

FAIRmessage([mfilename,' : done !']);

%==============================================================================
