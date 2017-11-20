%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
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
xn    = getNodalGrid(omega,m);

% setup interpolation scheme and image viewer
imgModel('reset','imgModel','linearInter');
viewImage('reset','viewImage','viewImage2D','colormap',gray(256),'axis','off');
fprintf('%20s : %s\n','viewImage',viewImage);
fprintf('%20s : %s\n','image model',imgModel);

% shortcuts for plotting stuff
Grid   = @(X)   plotGrid(X,omega,m,'spacing',8,'color','w');
dimstr = @(m)   sprintf('%s=[%s]',inputname(1),sprintf(' %d',m));
% Title  = @(s,t) title([s,sprintf(', t=%s',num2str(t))],'fontsize',20);


%% display initial image
% FAIRfigure(1,'figname',mfilename); clf; subplot(1,2,1); 
% vh = viewImage(imgModel(Tdata,omega,xc),omega,m); hold on; 
% gh = Grid(xc);
% title(sprintf('%s, %s','data',dimstr(m)),'fontsize',30);


yc = xn;
Tc = imgModel(Tdata,omega,center(yc,m));  

FAIRfigure(1,'figname',mfilename); clf; 
  subplot(1,2,1); %colordef(gcf,'black');
    vh = viewImage(Tc,omega,m); hold on; 
    gh = Grid(yc);
    title(sprintf('%s, %s','data',dimstr(m)),'fontsize',20);
  subplot(1,2,2); 
    vh = viewImage(Tc,omega,m); 
    th = title('model','fontsize',20);
FAIRpause(1);



dx = 0.3*(omega(2:2:end)-omega(1:2:end))
c  = (omega(2:2:end)+omega(1:2:end))'/2
S = @(t) diag([1-t/2,1-t/6]);
scale = @(t) reshape([S(t),(eye(2)-S(t))*c]',[],1);

zn = reshape(xn,[],2);
zn = [
  (zn(:,1)-omega(1))/(omega(2)-omega(1)),...
  (zn(:,2)-omega(3))/(omega(4)-omega(3))
  ];

zn = [ 
  0.5+(1-zn(:,2)).*cos(pi*(1-zn(:,1)))/2,...
  (1-zn(:,2)).*sin(pi*(1-zn(:,1)))
  ];

zn = [
  zn(:,1)*(omega(2)-omega(1))+omega(1)
  zn(:,2)*(omega(4)-omega(3))+omega(3)
  ];

cm = @(t) (1-t)*xn(:) + t*zn(:);
t  = sin(4*linspace(0,1,101)*pi);

trafos(1).name = 'translation2D-x';
trafos(1).map  = @(t) translation2D([dx(1)*t;0],xn);
trafos(2).name = 'translation2D-y';
trafos(2).map  = @(t) translation2D([0;dx(2)*t],xn);
trafos(3).name = 'translation2D-xy';
trafos(3).map  = @(t) translation2D([dx(1)*t;-dx(2)*t],xn);
trafos(4).name = 'rigid2D';
trafos(4).map  = @(t) rigid2D([t;(eye(2)-[ cos(t),-sin(t);sin(t),cos(t)])*c],xn);
trafos(5).name = 'scale';
trafos(5).map  = @(t) affine2D(scale(t),xn);
trafos(6).name = 'conformal';
trafos(6).map  = @(t) cm(t);
  

%=======================================================================================
% the loop over models and time 
%=======================================================================================

for trf = 1:length(trafos),
  set(th,'String',trafos(trf).name)
  for j=1:length(t),
    yc = trafos(trf).map(t(j));        % compute the transformed points
    Tc = imgModel(Tdata,omega,center(yc,m));   % compute the transformed image
    
    subplot(1,2,1);
    delete(gh);
    gh = Grid(yc);
    axis(omega);
    
    subplot(1,2,2);
    set(vh,'cdata',reshape(Tc,m)');
    FAIRpause(1/2000)
    
    drawnow; fprintf('.'); if rem(j,50) == 0, fprintf('\n'); end;
  end;
  fprintf('\n')
end;
 
%=======================================================================================
