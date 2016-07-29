% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% function coefT = getSplineCoefficients(dataT,varargin);
%
% computes regularized spline coefficients for interpolation, see FAIR Section 3.6.
% for 1D use coefT = getSplineCoefficients(dataT,'dim',1)
% version 2015/05/20
         
function coefT = getSplineCoefficients(dataT,varargin)

if nargin == 0, 
  runMinimalExample;
  return;
end;

% parameter and defaults
theta       = 0;
regularizer = 'none';
m           = size(dataT);
dim         = length(size(dataT));
out         = 1;

for k=1:2:length(varargin), % overwrite default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

if out, % report
  fprintf('compute %dD spline coefficients, m=[',dim);
  fprintf(' %d',m);
  fprintf('], regularizer=[%s], theta=[%s]\n',regularizer,num2str(theta));
end;

if isempty(dataT), coefT = []; return; end;
if ~strcmp(regularizer,'truncated'),
  switch dim,
    case 1,
      [M,P] = regularizedB(length(dataT),regularizer,theta);
      coefT = M\(P'*reshape(dataT,[],1));
    case 2,
      [M,P] = regularizedB(m(1),regularizer,theta);
      coefT = M\(P*dataT);
      coefT = permute(coefT,[2,1]);
      [M,P] = regularizedB(m(2),regularizer,theta);
      coefT = permute(M\(P*coefT),[2,1]);
    case 3,
      [M,P] = regularizedB(m(1),regularizer,theta);
      coefT = reshape(M\(P*reshape(dataT,m(1),m(2)*m(3))),m);
      [M,P] = regularizedB(m(2),regularizer,theta);
      coefT = reshape(permute(coefT,[2,1,3]),m(2),m(1)*m(3));
      coefT = permute(reshape(M\(P*coefT),m(2),m(1),m(3)),[2,1,3]);
      [M,P] = regularizedB(m(3),regularizer,theta);
      coefT = reshape(permute(coefT,[3,1,2]),m(3),m(1)*m(2));
      coefT = permute(reshape(M\(P*coefT),m(3),m(1),m(2)),[2,3,1]);
  end;
else
  switch dim,
    case 1,
      M = regularizedB(length(dataT),regularizer,theta);
      coefT = M*dataT;
    case 2,
      M = regularizedB(m(1),regularizer,theta);
      coefT = M*dataT;
      coefT = permute(coefT,[2,1]);
      M = regularizedB(m(2),regularizer,theta);
      coefT = permute(M*coefT,[2,1]);
    case 3,
      M = regularizedB(m(1),regularizer,theta);
      coefT = reshape(M*reshape(dataT,m(1),m(2)*m(3)),m);
      M = regularizedB(m(2),regularizer,theta);
      coefT = reshape(permute(coefT,[2,1,3]),m(2),m(1)*m(3));
      coefT = permute(reshape(M*coefT,m(2),m(1),m(3)),[2,1,3]);
      M = regularizedB(m(3),regularizer,theta);
      coefT = reshape(permute(coefT,[3,1,2]),m(3),m(1)*m(2));
      coefT = permute(reshape(M*coefT,m(3),m(1),m(2)),[2,3,1]);
  end;
end;
function [M,P] = regularizedB(m,regularizer,theta)
M = []; P = spdiags(ones(m,1)*[1,4,1],[-1:1],m,m);
switch regularizer,
  case 'none',      M = P;                        P = 1;
  case 'identity',  M = P'*P+theta*speye(m,m);
  case 'gradient',  D = spdiags(ones(m,1)*[-1,1],[0,1],m-1,m);
    M = P'*P+theta*D'*D;
  case 'moments',   M = P'*P+theta*toeplitz([96,-54,0,6,zeros(1,m-4)]);

  case 'truncated',
    k = max([1,min([ceil((1-theta)*m),m])]);
    d = 1./(4+2*cos((1:m)'*pi/(m+1)));
    V = sqrt(2/(m+1))*sin((1:m)'*(1:m)*pi/(m+1));
    M = V*diag([d(1:k);zeros(m-k,1)])*V';
  otherwise, error(regularizer);
end;
return;

function runMinimalExample

help(mfilename);
fprintf('%s: minimal example\n',mfilename)

omega = [0,10,0,8];
Tdata = [1,2,3,4;1,2,3,4;4,4,4,4;5,5,5,5]; m = size(Tdata);
Xdata = getCellCenteredGrid(omega,m);
xc    = getCellCenteredGrid(omega+[-1 1 -1 1],5*m);
Tcoef = getSplineCoefficients(Tdata,'out',1,'regularizer','moments','theta',1e-1)
[Tc,dT] = splineInter(Tcoef,omega,xc);
DD = reshape([Xdata;Tdata(:)],[],3);
Dc = reshape([xc;Tc],[5*m,3]);
figure(1); clf;
subplot(1,2,1); surf(Dc(:,:,1),Dc(:,:,2),Dc(:,:,3)); title(mfilename); 
subplot(1,2,2); spy(dT);                             title('dT')

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
