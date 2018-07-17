%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: SSD versus translations, splineInter, HNSP level=8
% 
%==============================================================================

clear, close all, help(mfilename)

setup2DHNSPData; 
level = 8; m = ML{level}.m; 
imgModel('set','imgModel','splineInter'); 
[T,R] = imgModel('coefficients',ML{level}.T,ML{level}.R,omega);
xc    = getCellCenteredGrid(omega,m); 
Rc    = imgModel(R,omega,xc);

trafo('set','trafo','translation2D');
figure(1); clf;
[w1,w2] = ndgrid(0.2*linspace(-1,1,21),0.2*linspace(-1,1,21));
dc  = zeros(size(w1));
for j=1:numel(dc),
  yc = trafo([w1(j);w2(j)],xc);
  Tc = imgModel(T,omega,yc);
  dc(j) = SSD(Tc,Rc,omega,m);
  viewImage(Tc,omega,m); FAIRpause(1/100)
end;
figure(1); clf; surf(w1,w2,dc); hold on; grid off; contour(w1,w2,dc)
title(sprintf('translation, m=[%d,%d]',m)); view(-135,33);
%==============================================================================
