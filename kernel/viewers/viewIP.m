%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
% ##2
%==============================================================================
% 
% function varargout = viewIP(I,omega,m,varargin)
%
% visualizes the intensity projections of a 3D image
%
% Input:
%   I           discretized image
%   omega		describing the domain
%	m 			number of discretization points
%   varargin    optional parameters like {'axis','off'}
%
% Output:
%  ih			image handle
%  I			[I1,I2;I3',0]
%==============================================================================

function varargout = viewIP(I,omega,m,varargin)

if nargin==0
help(mfilename);
    runMinimalExample; 
	return;
end

mode = 'max';
% reshape and form the average images
I  = reshape(I,m);
switch mode,
  case 'average',
    I1 = squeeze(sum(I,3))/m(3);
    I2 = squeeze(sum(I,2))/m(2);
    I3 = squeeze(sum(I,1))/m(1);
  case 'max',
    I1 = squeeze(max(I,[],3));
    I2 = squeeze(max(I,[],2));
    I3 = squeeze(max(I,[],1));
  case 'min',
    I1 = squeeze(min(I,[],3));
    I2 = squeeze(min(I,[],2));
    I3 = squeeze(min(I,[],1));
  otherwise, eror('nyi');
end;

I   = [I1,I2;I3',255+zeros(m(3),m(3))]; cla;
ih  = imagesc(I); axis image; hold on
ph1 = plot(0.5+m(2)*[1,1],0.5+[0,m(1)+m(3)],'b-','linewidth',2);
ph2 = plot(0.5+[0,m(2)+m(3)],0.5+m(1)*[1,1],'b-','linewidth',2);
ih  = [ih;ph1;ph2];

% the following lines add some nice stuff to the code.
% if varargin = {'title','FAIR','xlabel','x'}
% the code evaluates "title('FAIR');xlabel('x');"
for k=1:2:length(varargin), 
  if not(isempty(varargin{k})),
    feval(varargin{k},varargin{k+1}); 
  end;
end;

if nargout == 1, 
  varargout = {ih};
elseif nargout ==2,
  varargout = {ih,I};
end;

%------------------------------------------------------------------------------

function runMinimalExample
load mice3D; 
viewIP(dataT,omega,m);
%==============================================================================

