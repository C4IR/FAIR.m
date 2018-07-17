%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: 2D interpolation and transformations
%
% - load data ('USfair.jpg')
% - setup  viewer          (viewImage2D), 
%          interpolator    (linearInter), 
%          transformation  (rotation2D)
% - rotate image
%
% a non-trivial example for interpolation, transformtion, and visualization
% runs a loop over transformation parameters (rotation angles), 
%   computes the transformed grid, transformed image, and visualizes these
%==============================================================================

clear, close all, help(mfilename)

fprintf('%s\n','load data')
dataT = double(imread('US.jpg'));
m     = floor(size(dataT)/4);
omega = [0,size(dataT,1),0,size(dataT,2)];
xc    = getCellCenteredGrid(omega,m);

fprintf('%s\n','setup image viewer')
viewImage('reset','viewImage','viewImage2D','colormap',gray(256));

fprintf('%s\n','setup interpolator')
imgModel('reset','imgModel','linearInter');
TC = imgModel('coefficients',dataT,[],omega,'regularizer','moments','theta',1e1);
fprintf('%s\n','setup transformation model')
center = (omega(2:2:end)-omega(1:2:end))'/2;
trafo('reset','trafo','rotation2D','c',center);

fprintf('%s\n','start rotation')

% rotation angle wc shortcut for plots
wc  = sin(linspace(0,2*pi,201));
str = @(wc) sprintf('wc=%s',num2str(wc));
for j=1:length(wc),              % run over all angles
  yc = trafo(wc(j),xc);          % compute the transformed grid
  Tc = imgModel(TC,omega,yc);    % interpolate T on the grid

  if j == 1,                     % visualize the results
    FAIRfigure(1,'figname',mfilename);
    vh = viewImage(Tc,omega,m);
    th = title(str(wc(j)),'fontsize',20);
    set(gca,'fontsize',20);
    FAIRpause;
  else
    set(vh,'cdata',reshape(Tc,m)')
    set(th,'string',str(wc(j)))
    FAIRpause(1/2000)
  end;
  drawnow
  fprintf('.');
  if rem(j,50) == 0, fprintf('\n'); end;
end;
fprintf('\n');
%==============================================================================
