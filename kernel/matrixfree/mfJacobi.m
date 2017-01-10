%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function[x,r] = mfJacobi(x,b,para,steps)
% 
% Jacobi smoother for multigrid, see mfvcycle
%==============================================================================

function[x,r] = mfJacobi(x,b,para,steps)

if nargin == 0;
  MGsolver;
  x = 'endOfMinimalExample';
  return;
end;

% n = length(b);
[r,D] = mfAy(x,para);
r     = b - r; 

for i=1:steps,
  ss = D\r;
  x = x + para.MGomega*ss;   %x = rmnspace(x,para.Z);
  r = b - mfAy(x,para);      %r = b - rmnspace(mfAu(x,para),para.Z); 
%     his(i) = norm(r)/norm(b);  
%    figure(2); plot(r); pause
end;

%  figure(1); clf; plot(his); 
%  mfilename, keyboard
%==============================================================================
