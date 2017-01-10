%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: SSD versus rotations, linearInter, HNSP level=4
% 
%==============================================================================

clear, close all, help(mfilename);

setup2DHNSPData; 
imgModel('set','imgModel','linearInter'); 
level = 4; omega = ML{level}.omega; m = ML{level}.m; 
[T,R] = imgModel('coefficients',ML{level}.T,ML{level}.R,omega);
xc    = getCellCenteredGrid(omega,m); 
Rc    = imgModel(R,omega,xc);

center = (omega(2:2:end)-omega(1:2:end))'/2;
trafo('set','trafo','rotation2D','c',center);

wc = pi/2*linspace(-1,1,101);  dc = zeros(size(wc));
figure(1); clf;
for j=1:length(wc),
  yc = trafo(wc(j),xc);
  Tc = imgModel(T,omega,yc);
  dc(j) = SSD(Tc,Rc,omega,m);
  viewImage(128+(Tc-Rc)/2,omega,m); drawnow; 
  FAIRpause(1/6)
end;
figure; clf; p1 = plot(wc,dc); 
%==============================================================================
