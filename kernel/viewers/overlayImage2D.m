%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function varargout = overlayImage2D(T,R,omega,m,varargin)
%
% A 2D image viewer for overlaying images R and T
%
% Input:
%   T           discretized image
%   R           discretized image
%   omega       describing the domain
%    m          number of discretization points
%   varargin    optional parameters like {'colorT','r'}
%
% Output:
%  i1,i2        image handles for T and R
%==============================================================================

function varargout = overlayImage2D(T,R,omega,m,varargin)

if nargin==0
    runMinimalExample; 
	return;
end

colorT = [1,0,0];
colorR = [0,1,1];
scale  = 1;
alphadata = 0.5;

for k=1:2:length(varargin), % overwrites defaults
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

h  = (omega(2:2:end)-omega(1:2:end))./m;
xi = @(i) (omega(2*i-1)+h(i)/2:h(i):omega(2*i)-h(i)/2)';

if scale == 1,
  T = 255/max(T(:))*T;
  R = 255/max(R(:))*R;
end;


if length(size(T)) == 2,
  T = permute(reshape(uint8(T*colorT),[m,3]),[2,1,3]);
  R = permute(reshape(uint8(R*colorR),[m,3]),[2,1,3]);
else
  11
end;

cla;
i1  = image(xi(1),xi(2),R); axis xy image; hold on;
i2  = image(xi(1),xi(2),T); axis xy image; hold off;
set(i2,'alphaData',alphadata);

if nargout == 1, varargout = {[i1,i2]};  end;

%------------------------------------------------------------------------------
function runMinimalExample

help(mfilename)
setup2DhandData;
overlayImage2D(dataT(:),dataR(:),omega,m);
%==============================================================================

