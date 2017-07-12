%==============================================================================
% This code is part of the Matlab-based toolbox 
% LagLDDDM - A Lagrangian Gauss--Newton--Krylov Solver for Mass- and 
%                        Intensity-Preserving Diffeomorphic Image Registration
% 
% For details and license info see 
% - https://github.com/C4IR/FAIR.m/tree/master/add-ons/LagLDDMM
%
% function [ D ] = idctn( A, varargin )
%
% N-D Discrete Inverse Cosine Transform
%
% Computes the N-dimensional idct of an array, by applying the
% one-dimnesional idct along all dimensions.
% 
% Input:
%
%   A        -    full N-dimensional array
%   varargin - (optional flag in which dimensions to apply the idct)
%
% Output:
%
%   D        -    discrete inverse cosine transform of A
%
%==============================================================================

function [ D ] = idctn( A, varargin )

% size if the input array
m = size(A);
dim = length(m);

% default parameters
dimFlag = ones(1,dim);  % apply idct in all dimensions by default

% overwrites default parameter
for k=1:2:length(varargin),       
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

m = size(A);
D = A;
md = m;
P = circshift(1:dim,[0 -1]);

for d=1:dim
    if dimFlag(d)
        D = reshape(D,md(1),prod(md(2:end)));
        D = idct(D);
        D = reshape(D,md);
    end
    md = circshift(md,[0 -1]);
    D = permute(D,P);
end


