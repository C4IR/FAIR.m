% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% function [Tc,dT] = linearInterMatlab(T,omega,x,varargin)
%
% MATLAB's linear interpolator for the data T given on a cell-centered grid evaluated at x
% uses finite difference approximation for derivatives
% see Exercise 3.4
% version 2015/05/20

function [Tc,dT] = linearInterMatlab(T,omega,x,varargin)
         
Tc = mfilename('fullpath'); dT = []; 

if nargin == 0, 
  runMinimalExample;
  return;
elseif nargin == 1 && isempty(T),
  return;
end;
      
% if nargin == 0, return filename
if nargin == 0, Tc = mfilename('fullpath'); dT = []; return; end;


%  warning('not supported in MATLAB 2012')
dim = length(omega)/2;
n = numel(x)/dim;
Tc = zeros(n,1);
dT = sparse(n,dim*n);



% flag for computing the derivative
doDerivative = (nargout>1);
matrixFree   = 0;
for k=1:2:length(varargin), % overwrite default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

% get data size m, cell size h, dimension dim, and number n of interpolation points
dim = length(omega)/2;
m   = size(T);         if dim == 1, m = numel(T); end;
h   = (omega(2:2:end)-omega(1:2:end))./m; 
n   = length(x)/dim;    
x   = reshape(x,n,dim);

xi = @(i) (omega(2*i-1)+h(i)/2:h(i):omega(2*i)-h(i)/2)';
switch dim
    case 1
        xc = getCellCenteredGrid(omega,m);
        Tc = reshape(interp1(xc,T,x,'linear'),[],1);         Tc(isnan(Tc)) = 0;
        if ~doDerivative, return; end;
        h = (omega(2:2:end)-omega(1:2:end))./(length(T)*10);
        Tp = interp1(xc,T,x+h(1)/2,'linear');  Tp(isnan(Tp)) = 0;
        Tm = interp1(xc,T,x-h(1)/2,'linear');  Tm(isnan(Tm)) = 0;
        dT = (Tp-Tm)/h(1);
    case 2
        [X1,X2] = meshgrid(xi(1),xi(2));
        Tfctn = @(x1,x2) reshape(interp2(X1,X2,permute(T,[2,1]),x1,x2,'linear'),[],1);
        Tc = Tfctn(x(1:n),x(n+1:end));        Tc(isnan(Tc)) = 0; 
        if ~doDerivative, return; end;
        Tp = Tfctn(x(1:n)+h(1)/2,x(n+1:end)); Tp(isnan(Tp)) = 0;
        Tm = Tfctn(x(1:n)-h(1)/2,x(n+1:end)); Tm(isnan(Tm)) = 0;
        dT(:,1) = (Tp-Tm)/h(1);
        Tp = Tfctn(x(1:n),x(n+1:end)+h(2)/2); Tp(isnan(Tp)) = 0;
        Tm = Tfctn(x(1:n),x(n+1:end)-h(2)/2); Tm(isnan(Tm)) = 0;
        dT(:,2) = (Tp-Tm)/h(2);
    case 3
        [X1,X2,X3] = meshgrid(xi(1),xi(2),xi(3));
        Tfctn = @(x1,x2,x3) reshape(...
          interp3(X1,X2,X3,double(permute(T,[2,1,3])),x1,x2,x3,'linear'),[],1);
        Tc = Tfctn(x(1:n),x(n+1:2*n),x(2*n+1:end));         Tc(isnan(Tc)) = 0;
        if ~doDerivative, return; end;
        Tp = Tfctn(x(1:n)+h(1)/2,x(n+1:2*n),x(2*n+1:end));  Tp(isnan(Tp)) = 0;
        Tm = Tfctn(x(1:n)-h(1)/2,x(n+1:2*n),x(2*n+1:end));  Tm(isnan(Tm)) = 0;
        dT(:,1) = (Tp-Tm)/h(1);
        Tp = Tfctn(x(1:n),x(n+1:2*n)+h(2)/2,x(2*n+1:end));  Tp(isnan(Tp)) = 0;
        Tm = Tfctn(x(1:n),x(n+1:2*n)-h(2)/2,x(2*n+1:end));  Tm(isnan(Tm)) = 0;
        dT(:,2) = (Tp-Tm)/h(2);
        Tp = Tfctn(x(1:n),x(n+1:2*n),x(2*n+1:end)+h(3)/2);  Tp(isnan(Tp)) = 0;
        Tm = Tfctn(x(1:n),x(n+1:2*n),x(2*n+1:end)-h(3)/2);  Tm(isnan(Tm)) = 0;
        dT(:,3) = (Tp-Tm)/h(3);
end
if doDerivative && not(matrixFree)
    dT = spdiags(dT,n*(0:(dim-1)),n,dim*n);
end

function runMinimalExample
  help(mfilename);
  fprintf('%s: minimal examples\n',mfilename)

  % 1D example
  omega = [0,10];
  Tdata = [0,1,4,1,0]; 
  Tcoef = Tdata;
  m     = length(Tdata);
  Xdata = getCellCenteredGrid(omega,m);
  xc    = linspace(-1,11,101);
  [T0,dT0] = feval(mfilename,Tcoef,omega,xc);

  figure(1); clf;
  subplot(1,2,1); plot(xc,T0,'b-',Xdata,Tdata,'ro'); 
  title(sprintf('%s %d-dim',mfilename,1));
  subplot(1,2,2); spy(dT0);                     
  title('dT')

  % 2D example
  omega = [0,10,0,8];
  Tdata = [1,2,3,4;1,2,3,4;4,4,4,4]; m = size(Tdata);
  Tcoef = Tdata;
  Xdata    = getCellCenteredGrid(omega,m);
  xc    = getCellCenteredGrid(omega+[-1 1 -1 1],5*m);
  [Tc,dT] = feval(mfilename,Tcoef,omega,xc);
  DD = reshape([Xdata;Tdata(:)],[],3);
  Dc = reshape([xc;Tc],[5*m,3]);

  figure(2); clf;
  subplot(1,2,1);  surf(Dc(:,:,1),Dc(:,:,2),Dc(:,:,3));  hold on;
  plot3(DD(:,1),DD(:,2),DD(:,3),'r.','markersize',40); hold off;
  title(sprintf('%s %d-dim',mfilename,2));
  subplot(1,2,2); spy(dT);                     
  title('dT')

  % 3D example
  omega = [0,1,0,2,0,1]; m = [13,16,7];
  Xdata    = getCellCenteredGrid(omega,m);
  Y     = reshape(Xdata,[m,3]);
  Tdata = (Y(:,:,:,1)-0.5).^2 + (Y(:,:,:,2)-0.75).^2 + (Y(:,:,:,3)-0.5).^2 <= 0.15;
  Tcoef = double(reshape(Tdata,m));
  xc    = getCellCenteredGrid(omega,4*m);
  [Tc,dT] = feval(mfilename,Tcoef,omega,xc);

  figure(3); clf;
  subplot(1,2,1); imgmontage(Tc,omega,4*m);
  title(sprintf('%s %d-dim',mfilename,3));
  subplot(1,2,2); spy(dT);                 
  title('dT')

  fctn = @(xc) feval(mfilename,Tcoef,omega,xc);
  xc   = xc + rand(size(xc));
  checkDerivative(fctn,xc)
  
  %{ 
	=======================================================================================
	FAIR: Flexible Algorithms for Image Registration, Version 2011
	Copyright (c): Jan Modersitzki
	Maria-Goeppert-Str. 1a, D-23562 Luebeck, Germany
	Email: jan.modersitzki@mic.uni-luebeck.de
	URL:   http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
	=======================================================================================
	No part of this code may be reproduced, stored in a retrieval system,
	translated, transcribed, transmitted, or distributed in any form
	or by any means, means, manual, electric, electronic, electro-magnetic,
	mechanical, chemical, optical, photocopying, recording, or otherwise,
	without the prior explicit written permission of the authors or their
	designated proxies. In no event shall the above copyright notice be
	removed or altered in any way.

	This code is provided "as is", without any warranty of any kind, either
	expressed or implied, including but not limited to, any implied warranty
	of merchantibility or fitness for any purpose. In no event will any party
	who distributed the code be liable for damages or for any claim(s) by
	any other party, including but not limited to, any lost profits, lost
	monies, lost data or data rendered inaccurate, losses sustained by
	third parties, or any other special, incidental or consequential damages
	arrising out of the use or inability to use the program, even if the
	possibility of such damages has been advised against. The entire risk
	as to the quality, the performace, and the fitness of the program for any
	particular purpose lies with the party using the code.
	=======================================================================================
	Any use of this code constitutes acceptance of the terms of the above statements
	=======================================================================================
%}