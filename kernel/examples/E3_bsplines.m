%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: plots a B-spline
%
%==============================================================================

FAIRfigure(1,'color','w'); clf
tt = linspace(-3,11,1001);
p1 = plot(tt,spline1D(0,tt)); hold on; axis([-3,12,0,1])
set(p1,'linewidth',3,'color','k')
p2 = plot(tt,spline1D(2,tt),'--');
p3 = plot(tt,spline1D(7,tt),'--');
set([p2;p3],'linestyle','--','linewidth',1.5,'color','k'); axis off;
%==============================================================================
