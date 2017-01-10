%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function [t,Yt,LSiter] = Armijo(objFctn,Yc,dY,Jc,dJ,varargin)
%
% Armijo Line Search
%
% Armijo linesearch, returns t in [2^(-k), k=0:LSmaxIter], such that
%
%      objFctn( Yc + t*dY ) <= Jc + t*LSreduction*(dJ*dY) 
%
% (success!) or t=0 (no success)
%
% Input:
%   objFctn    function handle to the objective function
%   Yc         current vlue of Y
%   dY         search direction
%   Jc         current function value
%   dJ         current gradient 
%  varargin    optional parameters, see below
%
% Output:
%  t          steplength
%  Yt         new iterate
%  LSiter     number of steps performed
%
% see, e.g., 
%
%  @Book{NocWri1999,
%      author = {J. Nocedal and S. J. Wright},
%       title = {Numerical optimization},
%        year = {1999},
%   publisher = {Springer},
%     address = {New York},
%  }
%  
% see also GaussNewton and optPara.
%==============================================================================

function [t,Yt,LSiter,LS] = Armijo(objFctn,Yc,dY,Jc,dJ,varargin)

if nargin == 0, % help and minimal example
    help(mfilename);   
    GaussNewton;  
    t = 'endOfMinimalExample';
    return;
end;

LSmaxIter   = 10;           % max number of trials
LSreduction = 1e-4;         % slope of line
t = 1;                      % initial step
for k=1:2:length(varargin), % overwrites default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

L = [];
descent =   dJ * dY;
for LSiter =1:LSmaxIter,
  Yt = Yc + t*dY;       % compute test value Yt
  Jt = objFctn(Yt);      % evalute objective function
  LS = (Jt<Jc + t*LSreduction*descent); % compare
  if LS, break; end;    % success, return
  L = [L;t,Jt];
  t = t/2;          % reduce t
  
end;
if LS, return; end;      % we are fine
warning('Line Search failed - break\n');
t = 0; Yt = Yc;          % take no action

%==============================================================================
