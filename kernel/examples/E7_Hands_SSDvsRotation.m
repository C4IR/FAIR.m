%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: various distances versus rotation angle
%
%   - data                 xrays of hands, Omega=(0,20)x(0,25), level=7, m=[128,128]
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             plain L2-norm squared
%   - transformation       rotation2D
% see also E7_Hands_distance_rotation_ext
%==============================================================================

clear, close all, help(mfilename);

% load data, set viewer, interpolator, transformation, distance
setup2DhandData
viewImage('reset','viewImage','viewImage2D','colormap',bone(256),'axis','off');
imgModel('reset','imgModel','linearInter');
distance('reset','distance','SSD');
center = (omega(2:2:end)-omega(1:2:end))'/2;
trafo('reset','trafo','rotation2D','c',center);

% discretize T and R, make shortcut do difference image
xc = getCellCenteredGrid(omega,m);
Tc = imgModel(dataT,omega,xc);
Rc = imgModel(dataR,omega,xc);
Dc = @(Tc) 128+(Tc-Rc)/2;
hd = prod((omega(2:2:end)-omega(1:2:end))./m);

% discretize time and allocate memory
tt = linspace(0,2*pi,101);
D  = zeros(size(tt));
for k=1:length(tt),
  yc   = trafo(tt(k),xc);         % transform the grid X -> Y
  Tc   = imgModel(dataT,omega,yc);     % compute transformed image
  D(k) = hd*norm(Tc-Rc)^2;
  
  % visualize results
  if k == 1,    
    FAIRfigure(1,'figname',mfilename); 
    subplot(2,2,1); viewImage(Rc,omega,m);
    title('reference','fontsize',30)
    subplot(2,2,2); th = viewImage(Tc,omega,m);
    title('template','fontsize',30)
    subplot(2,2,3); dh = viewImage(Dc(Tc),omega,m);
    title('difference','fontsize',30)
    ax = [tt(1),tt(end),0,4e6];
    subplot(2,2,4); 
    plot(ax([1,1,2,2,1]),ax([3,4,4,3,3]),'w-','linewidth',2); hold on;
    axis(ax)
    plot(tt(1:k),D(1:k),'w-');
    ph = plot(tt(k),D(k),'r.','markersize',20);
    title('SSD versus rotation angle','fontsize',30)
    drawnow;
    FAIRpause;
  else
    set(th,'cdata',reshape(Tc,m)')
    set(dh,'cdata',reshape(Dc(Tc),m)')
    subplot(2,2,4); 
    plot(tt(1:k),D(1:k),'k-');
    set(ph,'visible','off');
    ph = plot(tt(k),D(k),'r.','markersize',20);
    drawnow;
    FAIRpause(1/2000)
  end;
end;

%==============================================================================
