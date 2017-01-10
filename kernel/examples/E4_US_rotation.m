%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR:  rotating an US image
%
%   - load US data
%   - set image viewer and interpolator and trafo=rigid2D
%   - run a loop over different angles
%==============================================================================

clear, close all, help(mfilename)

%% load data
Tdata = double(imread('US.jpg'));
omega = [0,size(Tdata,1),0,size(Tdata,2)];
m     = [192,128] ;
xc    = getCellCenteredGrid(omega,m);

% setup interpolation scheme and image viewer
imgModel('reset','imgModel','linearInter');
viewImage('reset','viewImage','viewImage2D','colormap',gray(256),'axis','off');
fprintf('%20s : %s\n','viewImage',viewImage);
fprintf('%20s : %s\n','image model',imgModel);

% shortcuts for plotting stuff
Grid   = @(X)   plotGrid(X,omega,m,'spacing',8,'color','w');
dimstr = @(m)   sprintf('%s=[%s]',inputname(1),sprintf(' %d',m));
Title  = @(s,t) title([s,sprintf(', t=%s',num2str(t))],'fontsize',30);


%% display initial image
FAIRfigure(1,'figname',mfilename); clf; subplot(1,2,1); 
viewImage(imgModel(Tdata,omega,xc),omega,m); hold on; gh = Grid(xc);
title(sprintf('%s, %s','data',dimstr(m)),'fontsize',30);

%% initialize transformation
trafo('set','trafo','rigid2D');
S = @(t) t/6;

% y(t) = R(t)*x+(I-R(t))*c, 
%
% R(t) = [cos(t) -sin(t) ], c = [omega(1)]/2
%        [sin(w)  cos(t) ]      [omega(2)]
%
% or y(t) = rigid2D(w,x), where w = [t;(I-R(t))omega'/2].
c  = (omega(2:2:end)+omega(1:2:end))'/2
wc = @(t) [t;(eye(2)-[ cos(t),-sin(t);sin(t),cos(t)])*c];
  
%=======================================================================================
% the loop over time t
%=======================================================================================
t = 0.35*pi*sin(linspace(0,2.5*pi,101));

for j=1:length(t),
  yc = trafo(wc(t(j)),xc);      % compute the transformed points
  Tc = imgModel(Tdata,omega,yc);   % compute the transformed image

  if j == 1,              % initialize visualization
    figure(1); subplot(1,2,1); %colordef(gcf,'black');
    set(gh,'visible','off'); gh = Grid(yc);
    title(sprintf('%s, %s','data',dimstr(m)),'fontsize',30);
    subplot(1,2,2); vh = viewImage(Tc,omega,m); Title(trafo,t(j));
    FAIRpause;
  else,                   % continue plots
    subplot(1,2,1); set(gh,'visible','off');        gh = Grid(yc);
    axis(omega)
    subplot(1,2,2); set(vh,'cdata',reshape(Tc,m)'); Title(trafo,t(j));
    FAIRpause(1/2000)
  end;
  drawnow; fprintf('.'); if rem(j,50) == 0, fprintf('\n'); end;
end;
fprintf('\n')
%=======================================================================================
