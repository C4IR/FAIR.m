%==============================================================================
% This code is part of the Matlab-based toolbox 
% LagLDDDM - A Lagrangian Gauss--Newton--Krylov Solver for Mass- and 
%                        Intensity-Preserving Diffeomorphic Image Registration
% 
% For details and license info see 
% - https://github.com/C4IR/FAIR.m/tree/master/add-ons/LagLDDMM
%
% function v0 = getVelocityStartingGuess(omega,m,nt)
%
% returns starting guess for velocity, that is, zero vector of appropiate
% size, which depends on regularizer, spatial and temporal discretization
% size.
%
% Input:
%
%  omega   - spatial domain
%      m   - number of cells for velocity
%     nt   - number of time points
%
% Output:
%
%     v0   - vector of zeros. size depends on grid for regularizer and nt
%
% =========================================================================
function v0 = getVelocityStartingGuess(omega,m,nt)
if not(exist('nt','var')) || isempty(nt)
    nt = regularizer('get','nt');% look up if regularizer is configured
    if isempty(nt),  nt = 0; end;
end
dim = numel(omega)/2;
switch regularizer('get','grid')
    case 'nodal'
        v0 = zeros(dim*prod(m+1),1);
    case 'staggered'
        v0 = zeros(sum(prod(ones(dim,1)*m+eye(dim),2)),1);
    otherwise % assume cell-centered
        v0 = zeros(dim*prod(m),1);
end
v0 = repmat(v0,nt+1,1);