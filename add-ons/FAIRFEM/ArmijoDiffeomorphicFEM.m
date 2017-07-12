%==============================================================================
% (c) Lars Ruthotto 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/lars-ruthotto.html
%
% function [t,Yt,LSiter,LS] = ArmijoBacktrackFEM(objFctn,Yc,dY,Jc,dJ,varargin)
%
% Armijo linesearch with additional volume control to guarantee a
% diffeomorphic update of the current iterate. 
%
% That is, the Armijo condition  of sufficient descent
%                        objFctn( Yc + t*dY ) <= Jc + t*LSreduction*(dJ*dY)
% is accompanied by the condition
%                        min(vol( Yc + t*dY )) > 0 
% which garantees that the nodal grid Yc + t*dY is diffeomorphic.
%
% if min(vol( Yc + t*dY ))>0 && objFctn( Yc + t*dY ) <= Jc + t*LSreduction*(dJ*dY), 
%   success!
% endIf
% t=2^-[0:10], else: t=0, no success
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
%  t      steplength
%  Yt      new iterate
%  LSiter    number of steps performed
%  LS      flag for success
%
% see, e.g., 
%  @Book{NocWri1999,
%      author = {J. Nocedal and S. J. Wright},
%       title = {Numerical optimization},
%        year = {1999},
%   publisher = {Springer},
%     address = {New York},
%  }
% and
%  @article{2011-BMR,
%	Author = {Burger M., Modersitzki J., Ruthotto L. },
%	Publisher = {University of Muenster},
%	Title = {A hyperelastic regularization energy for image registration},
%	Year = {2011}
%  }
%
% see also ArmijoBacktrack.m (version for nodal grids)
%==============================================================================

function [t,Yt,LSiter,LS] = ArmijoDiffeomorphicFEM(objFctn,Yc,dY,Jc,dJ,varargin)

if nargin == 0,
  help(mfilename)
  return;
end;

LSMaxIter   = 10;           % max number of trials
LSreduction = 1e-4;         % slope of line
para        = [];

for k=1:2:length(varargin), % overwrites default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

t = 1; descent =   dJ * dY; 
LS = 0; DIFFEOMORPHIC = 1;
for LSiter =1:LSMaxIter,
  Yt = Yc + t*dY; 			       % compute test value Yt
  V =  volTetraGrid(para.Mesh,Yt,'matrixFree',1);
  DIFFEOMORPHIC = (min(V(:))>0);    % check if update is diffeomorphic
  if DIFFEOMORPHIC,
    Jt = objFctn(Yt);              % evalute objective function
    LS = (Jt<Jc + t*LSreduction*descent); % compare
    if LS, break; end;             % success, return
  end
  t = t/2;                         % reduce t
end;
if LS, return; end;                % we are fine
if not(DIFFEOMORPHIC), 
    fprintf(['Line Search failed (No diffeomorphic update could be found)'...
        '- norm(dY)=%1.3e - break \n'],norm(dY));
elseif not(LS)
    fprintf(['Line Search failed (No sufficient descent found) '...
        '- norm(dY)=%1.3e - break\n'],norm(dY));
end
t = 0; Yt = Yc;        % take no action
