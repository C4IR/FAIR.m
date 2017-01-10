%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function [y,D] = mfAy(x,para)
%
% MatrixFree y =(M + alpha*h*B'*B)*y for elastic, diffusion, curvature
% where M is supposed to be diagonal and B is given via mfBu.m;
% the matrices M and B are represented by the struct para
% returns also D = diag(alpha*h*B'*B) 
%==============================================================================

function [y,D] = mfAy(x,para)

if nargin == 0,
  help(mfilename);
  E9_Hands_MLIR_SSD_mfElas;    
  y = 'endOfMinimalExample';
  return;
end;

omega = para.omega; 
m     = para.m; 

if isfield(para,'d2D') && isfield(para.d2D,'M')
  if isnumeric(para.d2D.M),
    Mx    = para.d2D.M.*x;
  elseif isa(para.d2D.M, 'function_handle')
      Mx = para.d2D.M(x);
  else
    error('M needs to be specified')
  end;
end;
y = Mx + para.d2S.d2S(x,omega,m);
if nargout<2, return; end;
if isfield(para,'d2D') && isfield(para.d2D,'M'),
    M = para.d2D.M;
else
    error('M needs to be specified');
end
% only for preconditioning of Jacobi smoothing in multigrid
S = para.d2S.diag(omega,m);
n = numel(S);
D = spdiags(M+S,0,n,n);
%==============================================================================
