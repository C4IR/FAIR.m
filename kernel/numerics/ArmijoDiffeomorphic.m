%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
% ##2
%==============================================================================
%
% function [t,Yt,LSiter] = ArmijoDiffeomorphic(objFctn,Yc,dY,Jc,dJ,varargin)
%
% Armijo linesearch with backtracking on the volume to guarantee diffeomorphic 
% updates for hyperelastic registration,   
% returns t in [2^(-k), k=0:LSmaxIter], such that
%	
%      objFctn( Yc + t*dY ) <= Jc + t*LSreduction*(dJ*dY) ...
%   && min(vol( Yc + t*dY )) > 0 
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
%------------------------------------------------------------------------------
% See, e.g., 
%
%  @Book{NocWri1999,
%      author = {J. Nocedal and S. J. Wright},
%       title = {Numerical optimization},
%        year = {1999},
%   publisher = {Springer},
%     address = {New York},
%  }
%
% @article{BurgerEtAl2013,
%    author = {Burger, M and Modersitzki, J and Ruthotto, L},
%    title = {{A hyperelastic regularization energy for image registration}},
%    journal = {SIAM Journal on Scientific Computing},
%    year = {2013},
%    volume = {35},
%    number = {1},
%    pages = {B132--B148},
%    doi = {10.1137/110835955},
% }
%  
% see also Armijo.m and Ehyper_2Ddisc2C.m.
%==============================================================================


function [t,Yt,LSiter,LS] = ArmijoBacktrack(objFctn,Yc,dY,Jc,dJ,varargin)

if nargin == 0,
   help(mfilename)
   Ehyper_2Ddisc2C_MB;
   t='end of minimal example';
   return;
end;

LSmaxIter   = 10;                       % max number of trials
LSreduction = 1e-4;                     % slope of line
para        = [];

for k=1:2:length(varargin),             % overwrites default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

m = para.m;

t = 1; descent =   dJ * dY; 
LS = 0; DIFFEOMORPHIC = 1;
for LSiter =1:LSmaxIter,
  Yt = Yc + t*dY; 			            % compute test value Yt
  range = min(geometry(Yt,m,'V'));      % compute distortion of volume
  diffeomorphic = (range>0);            % check if update is diffeomorphic
  if diffeomorphic,
    Jt = objFctn(Yt);                   % evalute objective function
    LS = (Jt<Jc + t*LSreduction*descent); % compare
    if LS, break; end;                  % success, return
  end
  t = t/2;                              % reduce t
end;
if LS, return; end;                     % we are fine

if not(diffeomorphic), 
  fprintf('%s // %s\n',...
	'Line Search failed : no diffeomorphic update ',...
    sprintf('norm(dY) = %1.3e',norm(dY)));
elseif not(LS)
  fprintf('%s // %s\n',...
	'Line Search failed : no sufficient descent',...
    sprintf('norm(dY) = %1.3e',norm(dY)));
end

t = 0; Yt = Yc;        % take no action
%==============================================================================
