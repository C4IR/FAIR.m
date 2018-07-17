%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: various transformations of an US image
%
%   - load US data
%   - set image viewer and interpolator
%   - visualize data
%   - run a loop over different transformation models (each model is parameterized by time t): 
%       translation-x1
%       translation-x2
%       translation
%       rotation
%       scale,
%       non-linear
%       spline
%==============================================================================

clear, close all, help(mfilename);

%% load US data
fprintf('%s\n','setup data')
dataT = double(imread('US.jpg'));
omega = [0,size(dataT,1),0,size(dataT,2)];
m     = [384,256]/2 ;
xc     = getCellCenteredGrid(omega,m);
center = (omega(2:2:end)+omega(1:2:end))'/2

fprintf('%s\n','setup image viewer and interpolator')
viewImage('reset','viewImage','viewImage2D','colormap',gray(256),'axis','off');
imgModel('reset','imgModel','linearInter'); T = dataT;
T0 = imgModel(T,omega,xc);
fprintf('%20s : %s\n','viewImage',viewImage);
fprintf('%20s : %s\n','imgModel',imgModel);

% shortcuts for plots
Grid    = @(xc) plotGrid(xc,omega,m,'spacing',16,'color','w');
dimstr  = @(m)  sprintf('%s=[%s]',inputname(1),sprintf(' %d',m));
Title  = @(s,t) title([s,sprintf(', t=%s',num2str(t))],'fontsize',30);

%%
fprintf('%s\n','visualize interpolated data')
figure(1); clf; colordef(gcf,'black'); 
subplot(1,2,1); viewImage(T0,omega,m); hold on; gh = Grid(xc);
title(sprintf('%s, %s','interpolated data',dimstr(m)),'fontsize',30);

% setup transformation models
trafos = {
  'translation-x1'
  'translation-x2'
  'translation'
  'rotation'
  'scale'
  'non-linear'
  'spline'
};

% discretize time
tt = sin(linspace(0,2.5*pi,201));

% pick a particular transformation, 
theTrafo = 1:length(trafos)

%=========================================================================================
%% a loop over transformations
%=========================================================================================

for k=theTrafo,

  fprintf('----- %s:\n',trafos{k});
  
  % based on the chosen transformation:
  %  setup trafo and the parameter wc as function in t
  clear trafo
  switch trafos{k},
    case 'translation-x1',              % translation in x1 direction
      trafo('set','trafo','translation2D');
      wc = @(t) [t*diff(omega(1:2))/4;0]; 
    case 'translation-x2',              % translation in x2 direction
      trafo('set','trafo','translation2D');
      wc = @(t) [0;t*diff(omega(3:4))/4]; 
    case 'translation',                 % translation in both directions
      trafo('set','trafo','translation2D');
      wc = @(t) [-t*diff(omega(1:2))/4;t*diff(omega(3:4))/6]; 
    case 'rotation',                    % rotation about the domain center
      trafo('set','trafo','rigid2D');
      S = @(t) t/6;
      wc = @(t) [S(t);(eye(2)-[ cos(S(t)),-sin(S(t));sin(S(t)),cos(S(t))])*center];
    case 'scale',                        % scale about the domain center
      trafo('set','trafo','affine2D');
      center = (omega(2:2:end)-omega(1:2:end))'/2;
      S = @(t) diag([1-t/4,1-t/6]);
      wc = @(t) reshape([S(t),(eye(2)-S(t))*center]',[],1);
    case 'non-linear',                   % a non-linear transformation
      % map omega to (-1,1) x (-1,1) and use a conformal map
      z  = (reshape(xc,[],2) - ones(prod(m),1)*center')*diag(0.5./center);
      z2 = [z.^2*[1;-1],2*z(:,1).*z(:,2)];
      trafo = @(t,xc) xc - 0.5*t*reshape(z2*diag(2*center),[],1);
      wc  = @(t) t;
    case 'spline',                   % a non-linear transformation
      p = [5,4];
      trafo('set','trafo','splineTransformation2D','p',p,'omega',omega,'m',m);
      w1 = zeros(p); w2 = zeros(p);  w2(3,2) = 10;
      wc  = @(t) [w1(:);t*w2(:)];        
    otherwise, error('nyi')
  end;%switch

  
  % the loop over time
  for j=1:length(tt),
    yc = trafo(wc(tt(j)),xc);   % transform the grid
    Tc = imgModel(T,omega,yc);     % compute transformed image
  
    % visualize the results
    if j == 1,                    
      FAIRfigure(k,'figname',mfilename); clf; 
      subplot(1,2,1);      viewImage(T0,omega,m); hold on; gh = Grid(yc);
      title(sprintf('%s, %s','interpolated data',dimstr(m)),'fontsize',30,'color','k');      
      subplot(1,2,2); vh = viewImage(Tc,omega,m);
      title(trafos{k},'fontsize',30,'color','k'); 
      FAIRpause;
    else
      subplot(1,2,1); set(gh,'visible','off'); gh = Grid(yc); axis(omega)
      subplot(1,2,2); set(vh,'cdata',reshape(Tc,m)'); FAIRpause(1/2000)
    end;
    drawnow; fprintf('.'); if rem(j,50) == 0, fprintf('\n'); end;
  end;
  fprintf('\n')
  FAIRpause;
end; % loop over transformations
%=========================================================================================


