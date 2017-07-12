%==============================================================================
% This code is part of the Matlab-based toolbox 
% LagLDDDM - A Lagrangian Gauss--Newton--Krylov Solver for Mass- and 
%                        Intensity-Preserving Diffeomorphic Image Registration
% 
% For details and license info see 
% - https://github.com/C4IR/FAIR.m/tree/master/add-ons/LagLDDMM
%
% function [dy,iter,relres,resvec] = spectralPrecondPCG(rhs, H, maxIterCG, tolCG,varargin)
%
% PCG with spectral preconditioner for the Gauss-Newton system in (MP)LDDMM as described 
% in the paper
%
% @article{MangRuthotto2017,
%   Title = {A {L}agrangian {G}auss--{N}ewton--{K}rylov solver for mass- and intensity-preserving diffeomorphic image registration},
%   Year = {2017},
%   Journal = {SIAM Journal on Scientific Computing},
%   Author = {A. Mang, L. Ruthotto},
% }
%
% Input: 
%
%    rhs       - vector, negative gradient
%    H         - struct describing Hessian
%    maxIterCG - max number of PCG iterations 
%    tolCG     - tolerance for PCG
%
% Optional Input (via varargin):
%
%    prec      - function handle for preconditioner, default='spectral'
%    shift     - modify Hessian by adding shift*speye(...) to ensure SPD
%
% Output:
%
%    dy        - search direction
%    iter      - iteration count
%    relres    - relative residual of returned solution
%    resvec    - PCG iteration history
% 
function [dy,iter,relres,resvec] = spectralPrecondPCG(rhs, H, maxIterCG, tolCG,varargin)

prec  = 'spectral';
shift = regularizer('get','HessianShift'); % shift Hessian approximation
%             to avoid rank-deficiency (similar to Levenberg - Marquardt)
if isempty(shift), shift = 0; end
for k=1:2:length(varargin),     % overwrites default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

% set pointers to simplify notation
omega   = H.omega;
m       = H.m;
tspan   = H.tspan;
nt      = H.nt;
P       = H.d2D.P; % grid-to-grid operator
dr      = H.d2D.dr;
d2psi   = H.d2D.d2psi;

% construct preconditioner
if isa(prec, 'function_handle')
    Preconditioner = prec;
else
    Preconditioner = @(x) solveSpectral(omega,m,tspan,nt,shift,1,x);
end

% build matrix free Hessian
Afctn = @(x) P(dr'*(d2psi*(dr*P(x)))) + H.d2S.d2S(x,H.omega,H.m) + shift*x;

% run PCG
[dy,flag,relres,iter,resvec] = pcg(Afctn,rhs,tolCG,maxIterCG,Preconditioner);

if exist('flag','var')
    switch flag
        case 1
            fprintf('pcg iterated %d times without converging to tolerance %e. Returned iterate (number %d) has residual %e\n',...
                numel(resvec)-1,tolCG,iter,relres);
            iter = numel(resvec)-1;
        case 2
            fprintf('Preconditioner M was ill-conditioned.\n');
        case 3
            fprintf('pcg stagnated. (Two consecutive iterates were the same.)\n');
        case 4
            fprintf('One of the scalar quantities calculated during pcg became too small or too large to continue computing.\n');
        otherwise
%             fprintf('pcg success! %d iterations / relres= %1.2e / tolCG= %1.2e\n',iter,relres,tolCG);
    end
end
