%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function function yc = center(yc,m,varargin)
%
% transfers any yc (centered, staggered, nodal) to a centered grid
%
% Input:
%   yc      current grid points
%   m       number of discretization points
%
% Output:
%   yc      on cell-centered grid
%
% for 1D: yc = center(yc,m,'dim',1);
% uses nodal2center, stg2center
%==============================================================================

function yc = center(yc,m,varargin)

if nargin==0
    help(mfilename)
    runMinimalExample; 
    return;
end

dim = length(m);
for k=1:2:length(varargin),     % overwrites default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

if dim*prod(m) == numel(yc),
  % yc is already cell-centered
  yc = reshape(yc,[],1);
elseif dim*prod(m+1) == numel(yc),
  % yc is nodal 
  yc = nodal2center(yc,m);
elseif sum(prod(ones(dim,1)*m+eye(dim),2)) == numel(yc)
  % yc is staggered
  yc = stg2center(yc,m);
else
  error('don''t know how to deal this grid')
end;

%------------------------------------------------------------------------------

function runMinimalExample
fprintf('Example:\n')
yNodal  = (0:2:10)
yCenter = center(yNodal,length(yNodal)-1,'dim',1)'

%==============================================================================
