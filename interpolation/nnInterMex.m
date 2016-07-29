% (c) Fabian Gigengack 2011/08/12, see FAIR.2 and FAIRcopyright.m.
% http://www.uni-muenster.de/EIMI/
%
% function [Tc,dT] = nnInterMex(T,omega,x,varargin)
%
% This wrapper runs a CPP file which codes the same scheme as nnInter
% but is hopefully faster, see nnInter for details
% version 2015/05/20

function [Tc,dT] = nnInterMex(T,omega,x,varargin)

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
boundary     = 0; % boundary condition (0: zero padding; 1: replicate)

for k=1:2:length(varargin) % overwrite default parameter
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end

dim = length(omega)/2;
m   = size(T);         if dim == 1, m = numel(T); end;
n   = length(x)/dim;

% call CPP subroutine
try
    Tc = nnInterMexC(double(T(:)),omega,m,x(:),boundary==1);
catch err
    FAIRerror(err);
end
if doDerivative
    if matrixFree, dT = zeros(n,dim); else dT = sparse(n,dim*n); end
end

function runMinimalExample
help(mfilename);
fprintf('%s: minimal examples\n',mfilename)

% 1D example
omega = [0,10];
TD    = [0,1,4,1,0];
m     = length(TD);
XD    = getCellCenteredGrid(omega,m);
xc    = linspace(-1,11,101);
[T0,dT0] = feval(mfilename,TD,omega,xc);

figure(1);
subplot(1,2,1); plot(xc,T0,'b-',XD,TD,'ro');
title(sprintf('%s %d-dim',mfilename,1));
subplot(1,2,2); spy(dT0);
title('dT')

% 2D example
omega = [0,10,0,8];
TD    = [1,2,3,4;1,2,3,4;4,4,4,4]; m = size(TD);
XD    = getCellCenteredGrid(omega,m);
xc    = getCellCenteredGrid(omega+[-1 1 -1 1],5*m);
[Tc,dT] = feval(mfilename,TD,omega,xc);
DD = reshape([XD;TD(:)],[],3);
Dc = reshape([xc;Tc],[5*m,3]);

figure(2); clf;
subplot(1,2,1);  surf(Dc(:,:,1),Dc(:,:,2),Dc(:,:,3));  hold on;
plot3(DD(:,1),DD(:,2),DD(:,3),'r.','markersize',40); hold off;
title(sprintf('%s %d-dim',mfilename,2));
subplot(1,2,2); spy(dT);
title('dT')

% 3D example
omega = [0,1,0,2,0,1]; m = [13,16,7];
XD    = getCellCenteredGrid(omega,m);
Y     = reshape(XD,[m,3]);
TD    = (Y(:,:,:,1)-0.5).^2 + (Y(:,:,:,2)-0.75).^2 + (Y(:,:,:,3)-0.5).^2 <= 0.15;
TD    = reshape(TD,m);
xc    = getCellCenteredGrid(omega,4*m);
[Tc,dT] = feval(mfilename,TD,omega,xc);

figure(3); clf;
subplot(1,2,1); imgmontage(Tc,omega,4*m);
title(sprintf('%s %d-dim',mfilename,3));
subplot(1,2,2); spy(dT);
title('dT')

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