% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% function b = motherSpline(j,xi)
%
% shortcut for the mother spline function on the four non-trivial intervals (1,2,3,4) 
% and its derivatives (5,6,7,8), see FAIR, p.27 Eq. (3.7).
% version 2015/05/20

function b = motherSpline(j,xi)

if nargin == 0,
  help(mfilename);
  fprintf('%s: minimal example\n',mfilename)
  m = 101;  
  x = linspace(0,1,m);
  z = zeros(m,4);  
  y = zeros(m,4,2);  
  for j=1:4,
    z(:,j)   = x+j-3;
    y(:,j,1) = motherSpline(j,  z(:,j)); 
    y(:,j,2) = motherSpline(j+4,z(:,j)); 
  end;
  figure(1); clf; 
  plot(z(:),reshape(y(:,:,1),1,[]),'b-',...
    z(:),reshape(y(:,:,2),1,[]),'r-'); 
  title(mfilename);
  legend('motherspline','derivative')
  %plot(x(1,:),x(2:end,:)); title(mfilename);
  return;
end;

switch j,
  case 1,  b = (2+xi).*(2+xi).*(2+xi);  % MATLAB! b = (2+xi).^3;
  case 2,  b = -(3*xi+6).*xi.^2+4;      % -xi.^3 - 2*(xi+1).^3 +6*(xi+1);
  case 3,  b =  (3*xi-6).*xi.^2+4;      % xi.^3 + 2*(xi-1).^3 -6*(xi-1);
  case 4,  b = (2-xi).*(2-xi).*(2-xi);  % MATLAB! b = (2-xi).^3;
  case 5,  b =  3*(2+xi).^2;
  case 6,  b = -(9*xi+12).*xi;          % -3*xi.^2 - 6*(xi+1).^2+6;
  case 7,  b =  (9*xi-12).*xi;          %  3*xi.^2 + 6*(xi-1).^2-6;
  case 8,  b = -3*(2-xi).^2;
end;

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