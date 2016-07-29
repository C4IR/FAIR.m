% (c) Fabian Gigengack 2011/04/13 see FAIR.2 and FAIRcopyright.m.
% http://www.uni-muenster.de/EIMI/
%
% function [Tc,dT] = cubicInter(T,omega,x,varargin)
%
% Fast cubic interpolator with Catmull-Rom splines for the data T given on
% a cell-centered grid on omega, to be evaluated at x
% 
% See:
% Catmull, E., and Rom,  R. A class of local interpolating splines.
% In Computer Aided Geometric Design, R. E. Barnhill and R. F. Reisenfeld,
% Eds. Academic Press, New York, 1974, pp. 317-326.

function [Tc,dT] = cubicInter(T,omega,x,varargin)
         
% if nargin == 0, return filename
Tc = mfilename('fullpath'); dT = [];

if nargin == 0
    runMinimalExample;
    return
elseif nargin == 1 && isempty(T)
    return
end

% flag for computing the derivative
doDerivative = (nargout>1);
matrixFree   = 0;
for k=1:2:length(varargin), % overwrite default parameter
	eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end

% get data size m, cell size h, dimension d, and number n of interpolation points
d = length(omega)/2;
m = size(T); if d==1, m = numel(T); end
h = (omega(2:2:end)-omega(1:2:end))./m;
n = length(x)/d;
x = reshape(x,n,d);
% map x from [h/2,omega-h/2] -> [1,m],
for i=1:d, x(:,i) = (x(:,i)-omega(2*i-1))/h(i) + 0.5; end

vec = @(in) in(:);
Tc = zeros(n,1); dT = [];                   % initialize output
if doDerivative, dT = zeros(n,d);  end      % allocate memory in column format
Valid = @(j) (0<x(:,j) & x(:,j)<m(j)+1);    % determine indices of valid points

switch d
    case 1, valid = find(Valid(1));
    case 2, valid = find(Valid(1) & Valid(2));
    case 3, valid = find(Valid(1) & Valid(2) & Valid(3));
end

if isempty(valid),                        
    if doDerivative, dT = sparse(n,d*n); end % allocate memory incolumn format
    return
end

pad = 2;
P  = floor(x); x = x-P;                   % split x into integer/remainder
p  = @(j) P(valid,j); xi = @(j) x(valid,j);
i1 = 1; i2 = size(T,1)+2*pad; i3 = (size(T,1)+2*pad)*(size(T,2)+2*pad);

switch d
    case 1
        TP = zeros(m+2*pad,1);
        TP(pad+(1:m)) = reshape(T,m,1);
        clear T;
        p = pad + P(valid); xi = x(valid); % add the padding
        idx = @(i) p+(i-2);
        % cubic interpolation in x-dimension
        Tc(valid) = cint(TP(idx(1)), TP(idx(2)), TP(idx(3)), TP(idx(4)), vec(xi));
    case 2
        TP = zeros(m+2*pad);
        TP(pad+(1:m(1)),pad+(1:m(2))) = T;
        clear T;
        p   = (pad + p(1)) + i2*(pad + p(2) - 1);
        idx = @(i,j) p+(i-2)*i1+(j-2)*i2;
        % cubic interpolation in y-dimension
        u = zeros(4,numel(p));
        for i=1:4 % x-dimension
            u(i,:) = cint(TP(idx(i,1)), TP(idx(i,2)), TP(idx(i,3)), TP(idx(i,4)), vec(xi(2)));
        end
        % cubic interpolation in x-dimension
        Tc(valid) = cint(u(1,:)', u(2,:)', u(3,:)', u(4,:)', vec(xi(1)));
    case 3
        TP = zeros(m+2*pad);
        TP(pad+(1:m(1)),pad+(1:m(2)),pad+(1:m(3))) = T;
        clear T;
        p   = (pad + p(1)) + i2*(pad + p(2) - 1) + i3*(pad + p(3) -1);
        idx = @(i,j,k) p+(i-2)*i1+(j-2)*i2+(k-2)*i3;
        % cubic interpolation in z-dimension
        t = zeros(4,4,numel(p));
        for i=1:4 % x-dimension
            for j=1:4 % y-dimension
                t(i,j,:) = cint(TP(idx(i,j,1)), TP(idx(i,j,2)), ...
                                TP(idx(i,j,3)), TP(idx(i,j,4)), vec(xi(3)));
            end
        end
        % cubic interpolation in y-dimension
        u = zeros(4,numel(p));
        for i=1:4
            u(i,:) = cint(vec(t(i,1,:)), vec(t(i,2,:)), vec(t(i,3,:)), vec(t(i,4,:)), vec(xi(2)));
        end
        % cubic interpolation in x-dimension
        Tc(valid) = cint(u(1,:)', u(2,:)', u(3,:)', u(4,:)', vec(xi(1)));
end

if ~doDerivative, return; end

switch d
    case 1
        % compute and format the derivative
        dT(valid) = cint(TP(idx(1)), TP(idx(2)), TP(idx(3)), TP(idx(4)), vec(xi), true);
    case 2
        % derivative of cubic interpolation in x-direction
        dT(valid,1) = cint(u(1,:)', u(2,:)', u(3,:)', u(4,:)',  vec(xi(1)), true);
        % derivative of cubic interpolation in y-direction
        u = zeros(4,numel(p));
        for i=1:4
            u(i,:) = cint(TP(idx(i,1)), TP(idx(i,2)), TP(idx(i,3)), TP(idx(i,4)),  vec(xi(2)), true);
        end
        dT(valid,2) = cint(u(1,:)', u(2,:)', u(3,:)', u(4,:)', vec(xi(1)));
    case 3
        % derivative of cubic interpolation in x-direction
        dT(valid,1) = cint(u(1,:)', u(2,:)', u(3,:)', u(4,:)', vec(xi(1)), true);
        % derivative of cubic interpolation in y-direction
        u = zeros(4,numel(p));
        for i=1:4
            u(i,:) = cint(vec(t(i,1,:)), vec(t(i,2,:)), vec(t(i,3,:)), vec(t(i,4,:)), vec(xi(2)), true);
        end
        dT(valid,2) = cint(u(1,:)', u(2,:)', u(3,:)', u(4,:)', vec(xi(1)));
        % derivative of cubic interpolation in z-direction
        t = zeros(4,4,numel(p));
        for i=1:4
            for j=1:4
                t(i,j,:) = cint(TP(idx(i,j,1)), TP(idx(i,j,2)), ...
                                TP(idx(i,j,3)), TP(idx(i,j,4)), vec(xi(3)), true);
            end
        end
        u = zeros(4,numel(p));
        for i=1:4
            u(i,:) = cint(vec(t(i,1,:)), vec(t(i,2,:)), vec(t(i,3,:)), vec(t(i,4,:)), vec(xi(2)));
        end
        dT(valid,3) = cint(u(1,:)', u(2,:)', u(3,:)', u(4,:)', vec(xi(1)));
end
for i=1:d, dT(:,i) = dT(:,i)/h(i); end
if not(matrixFree)
    dT = spdiags(dT,n*(0:(d-1)),n,d*n);
end

function out = cint(p0, p1, p2, p3, v, doDerivative)

if nargin<6, doDerivative = false; end

v2 = v.^2; v3 = v.^3;

if ~doDerivative
    % Cubic interpolation function
    out =   (- .5 * v3 +       v2 - .5 * v    ) .* p0 ...
          + ( 1.5 * v3 - 2.5 * v2          + 1) .* p1 ...
          + (-1.5 * v3 + 2   * v2 + .5 * v    ) .* p2 ...
          + (  .5 * v3 -  .5 * v2             ) .* p3;
else
    % Derivative of cubic interpolation function
    out =   (-1.5 * v2 + 2 * v - .5) .* p0 ...
          + ( 4.5 * v2 - 5 * v     ) .* p1 ...
          + (-4.5 * v2 + 4 * v + .5) .* p2 ...
          + ( 1.5 * v2 -     v     ) .* p3;
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
  Tcoef = reshape(Tdata,m);
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