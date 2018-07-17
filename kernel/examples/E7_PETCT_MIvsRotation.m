%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: MI versus rotation
%
%   - data                 PETCT, Omega=(0,140)x(0,151), level=4:7, m=[128,128]
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             MI
%   - transformation       rotation2D
%==============================================================================

clear, close all, help(mfilename);

% setup data, interpolation, transformation, regularization
setup2DPETCTData; level = 6; omega = ML{level}.omega; m = ML{level}.m;

viewImage('reset',viewPara{:},'axis','off');
imgModel('reset','imgModel','splineInter','regularizer','moments','theta',1e0);
[T,R] = imgModel('coefficients',ML{level}.T,ML{level}.R,omega,'out',0);
distance('reset','distance','MI');
center = (omega(2:2:end)-omega(1:2:end))'/2;
trafo('reset','trafo','rotation2D','c',center);
fprintf('%20s : %s\n','viewImage',viewImage);
fprintf('%20s : %s\n','imgModel',imgModel);
fprintf('%20s : %s\n','distance',distance);
fprintf('%20s : %s\n','trafo',trafo);

xc = getCellCenteredGrid(omega,m);
Rc = imgModel(R,omega,xc);

% diffImage = @(Tc) viewImage(128+(Tc-Rc)/2,omega,m);
wc = linspace(-pi/2,pi/2,51);
Dc = zeros(size(wc));

% run the loop over all rotations
for j = 1:length(wc),
  yc = trafo(wc(j),xc);                 % compute transformed grid
  Tc = imgModel(T,omega,yc);               % compute transformed image
  [Dc(j),rc] = distance(Tc,Rc,omega,m); % compute distance
  
  % visualize
  if j == 1, th = [];
    FAIRfigure(1,'figname',mfilename); clf;
    subplot(1,3,1);  viewImage(Rc,omega,m);            th(1) = title('R');
    subplot(1,3,2);  vh = viewImage(Tc,omega,m);       th(2) = title('T(yc)');
    subplot(1,3,3);  ph = plot(wc(1),Dc(1),'r.','markersize',20);
    th(3) = title(sprintf('%s versus rotation',distance));
    axis([wc(1),wc(end),-inf,inf]); hold on;
    axis('auto y')
    set(th,'fontsize',30);
    FAIRpause;
  else
    set(vh,'cdata',reshape(Tc,m)')
    subplot(1,3,3); set(ph,'visible','off');
    plot(wc(1:j),Dc(1:j),'k-','linewidth',2);
    ph = plot(wc(j),Dc(j),'r.','markersize',20);
    drawnow, 
    FAIRpause(1/100)
  end;
  fprintf('.'); if ~rem(j,50) || j == length(wc), fprintf('\n'); end;
end;

%==============================================================================
