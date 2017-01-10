%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: US, MI versus rotation
%
%   - data                 US
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             MI
%   - transformation       rotation2D
%   - optimization         Gauss-Newton
% ===============================================================================

clear, close all, help(mfilename);

% load data
dataT = double(imread('US.jpg'));
dataR = 255 - dataT;
omega = [0,size(dataT,1),0,size(dataT,2)];
m     = [192,128]/2;
xc    = getCellCenteredGrid(omega,m);

% setup image viewer, interpolation, transformation, distance
viewImage('reset','viewImage','viewImage2D','colormap',gray(256));
imgModel('reset','imgModel','splineInter','regularizer','moments','theta',1e0);
% don't forget to compute the coefficients!
[T,R] = imgModel('coefficients',dataT,dataR,omega,'out',0);
distance('reset','distance','MI',...
  'tol',1e-7,'minT',0,'maxT',256,'nT',64,'minR',0,'maxR',256,'nR',64);

center = (omega(2:2:end)-omega(1:2:end))'/2;
trafo('reset','trafo','rotation2D','c',center);

fprintf('%20s : %s\n','viewImage',viewImage);
fprintf('%20s : %s\n','imgModel',imgModel);
fprintf('%20s : %s\n','trafo',trafo);
fprintf('%20s : %s\n','distance',distance);
distance('disp');

wc = linspace(-pi/2,pi/2,101);
Dc = zeros(size(wc));
Rc = imgModel(R,omega,xc);

% prepare plots  
nT  = distance('get','nT') + 4; % reformat and scale the joint density rho
nR  = distance('get','nR') + 4;
rho = @(rc) reshape((1+min(rc,0.001)),nT,nR);
  
% run the loop over all rotations
for j = 1:length(wc),
  yc = trafo(wc(j),xc);        % compute transformed points
  Tc = imgModel(T,omega,yc);    % compute transformed image
  [Dc(j),rc] = distance(Tc,Rc,omega,m); % compute distance and joint density

  if j == 1,    % initialize plots
    FAIRfigure(1,'figname',mfilename); clf;
    subplot(2,2,1); viewImage(Rc,omega,m);             th(1) = title('R');
    subplot(2,2,2); vh = viewImage(Tc,omega,m);        th(2) = title('T(y)');
    subplot(2,2,3); rh = imagesc(rho(rc)'); axis image; th(3) = title('\rho(T(yc),R)');
    subplot(2,2,4); ph = plot(wc(j),Dc(j),'r.','markersize',20);
    th(4) = title(sprintf('%s versus rotation',distance));
    set(th,'fontsize',30); hold on;
    ax = [wc(1),wc(end),50,60];  
    FAIRpause;
  else,         % update plot
    figure(1); set(vh,'cdata',reshape(Tc,m)'); set(rh,'cdata',rho(rc)');
    subplot(2,2,4); set(ph,'visible','off');
    plot(wc(1:j),Dc(1:j),'r-','linewidth',2);
    ph = plot(wc(j),Dc(j),'r.','markersize',20);
    drawnow, 
    FAIRpause(1/1000);
  end;
  fprintf('.');
  if ~rem(j,50) || j == length(wc), fprintf('\n'); end;
end;
%==============================================================================
