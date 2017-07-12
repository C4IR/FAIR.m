%==============================================================================
% ##1
%
% function [ D ] = dctn( A, varargin )
%
% N-D Discrete Cosine Transform
%
% Computes the N-dimensional dct of an array, by applying the
% one-dimnesional dct along all dimensions.
% 
% Input:
%
%   A        - full N-dimensional array
%   varargin - (optional flag in which dimensions to apply the dct)
%
% Output:
%
%   D        - discrete cosine transform of A
%
%==============================================================================

function [ D ] = dctn( A, varargin )

if nargin==0
    help(mfilename)
    return;
end

% size if the input array
m = size(A);
dim = length(m);

% default parameters
dimFlag = ones(1,dim);  % apply dct in all dimensions by default

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
        D = dct(D);
        D = reshape(D,md);
    end
    md = circshift(md,[0 -1]);
    D = permute(D,P);
end


