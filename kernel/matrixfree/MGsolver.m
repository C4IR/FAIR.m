%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
% ##2
%==============================================================================
%
% function Y = MGsolver(rhs,H);
%
% initializes a multigrid call (mfvcycle)  for H*Y=rhs
%==============================================================================

function Y = MGsolver(rhs,H,varargin)

if nargin == 0,
  help(mfilename);
  E9_Hands_MLIR_SSD_mfElas;
  Y = 'endOfMinimalExample';
  return;
end;

MGlevel      = log2(H.m(1))+1;%%(~strcmp(H.d2S.regularizer,'mfCurvature'));
MGcycle      = 1;
MGomega      = 0.5;   
MGsmoother   = 'mfJacobi';
MGpresmooth  = 3;
MGpostsmooth = 1;

for k=1:2:length(varargin), % overwrites default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

H.MGlevel      = MGlevel;     
H.MGcycle      = MGcycle;     
H.MGomega      = MGomega;     
H.MGsmoother   = MGsmoother;  
H.MGpresmooth  = MGpresmooth; 
H.MGpostsmooth = MGpostsmooth;

d2D         = H.d2D;
H.d2D       = [];


if isfield(d2D,'M'), % get approximation to d2D
  M     = d2D.M;
elseif isfield(d2D,'dr')
  if isscalar(d2D.d2psi)
     M           =  sum(d2D.dr .* d2D.dr,1) .* d2D.d2psi;
  else
      M           = diag(d2D.dr'*d2D.d2psi*d2D.dr);
  end
  M           = d2D.P(full(M(:)));
else
  M = 0;
end;

% bring M into right format
if all(size(M)==1),
    M = M * ones(length(rhs),1);
elseif (any(size(M)) == 1)
    M = full(M(:));
else
    error('M needs to be diagonal and stored as a vector!');
end
H.d2D.M = M;
    
u0 = zeros(size(rhs));
Y  = mfvcycle(H,u0,rhs,1e-12,H.MGlevel,(max(H.m)>32));

%==============================================================================
