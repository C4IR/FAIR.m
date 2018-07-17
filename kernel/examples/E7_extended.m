%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: various distances versus rotation angle
%
%   - data                 PETCT, Omega=(0,140)x(0,151), level=4:7, m=[128,128]
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             {'SSD','NCC','MI','NGF'}
%   - transformation       rotation2D
% see also E7_basic
%==============================================================================

clear, close all, help(mfilename);

setup2DPETCTData; level = 6; omega = ML{level}.omega; m = ML{level}.m;
viewImage('reset',viewPara{:},'axis','off');
imgModel('reset','imgModel','splineInter','regularizer','moments','theta',1e0);
[T,R] = imgModel('coefficients',ML{level}.T,ML{level}.R,omega,'out',0);
distance('reset','distance','SSD');
center = (omega(2:2:end)-omega(1:2:end))'/2;
trafo('reset','trafo','rotation2D','c',center);
fprintf('%20s : %s\n','viewImage',viewImage);
fprintf('%20s : %s\n','imgModel',imgModel);
fprintf('%20s : %s\n','distance',distance);
fprintf('%20s : %s\n','trafo',trafo);

xc = getCellCenteredGrid(omega,m);
Rc = imgModel(R,omega,xc);
wc = linspace(-pi/2,pi/2,51);
Dc = zeros(size(wc));

distances = {'SSD','NCC','MI','NGF'};
diffImage = @(Tc) viewImage(128+(Tc-Rc)/2,omega,m);

for k=1:length(distances),

  distance('reset','distance',distances{k});
  fprintf('---- %s\n\n',distance);
  
  Dc = zeros(size(wc));
  % adapt visualization functions
  switch distance,
    case 'SSD',
      resStr = '|T(yc)-R|';
      res = @(rc) 128+rc/2;
      plotRes = @(rc) viewImage(res(rc),omega,m);
      setRes  = @(ph,rc) set(ph,'cdata',reshape(res(rc),m)');
    case 'NCC',
      resStr = '|T(yc)-R|';
      res = @(rc) 128+rc/2;
      plotRes = @(rc) viewImage(res(rc),omega,m);
      setRes  = @(ph,rc) set(ph,'cdata',reshape(res(rc),m)');
    case 'MI',
      %setup default parameter
      distance('set','tol',1e-7,'minT',0,'maxT',256,'nT',60,'minR',0,'maxR',256,'nR',60);
      distance('disp');
      nT = distance('get','nT') + 4;
      nR = distance('get','nR') + 4;
      ax = [wc(1),wc(end),50,60];
      resStr = '\rho(|T(yc),R)';
      res = @(rc) reshape((1+min(rc,0.001)),nT,nR);
      plotRes = @(rc) imagesc(res(rc));
      setRes  = @(ph,rc) set(ph,'cdata',res(rc)');
    case 'NGF',
      %setup default parameter
      distance('set','edge',100);
      distance('disp');
      ax = [wc(1),wc(end),370,400];
      resStr = 'NGF(T(yc),R)';
      res = @(rc) reshape(rc,m);
      plotRes = @(rc) viewImage2Dsc(res(rc),omega,m);
      setRes  = @(ph,rc) set(ph,'cdata',res(rc)');
  end;

  % run the loop over all rotations
  for j = 1:length(wc),
    yc = trafo(wc(j),xc);
    Tc = imgModel(T,omega,yc);
    [dc,rc] = distance(Tc,Rc,omega,m);
    Dc(j) = dc;

    if j == 1,
      FAIRfigure(k,'figname',mfilename); clf;
      subplot(2,2,1);  viewImage(Rc,omega,m);          th(1) = title('R');
      subplot(2,2,2);  vh = viewImage(Tc,omega,m);     th(2) = title('T(yc)');
      subplot(2,2,3);  rh = plotRes(rc); axis image;   th(3) = title(resStr);
      subplot(2,2,4);  ph = plot(wc(1),Dc(1),'r.','markersize',20);
      th(4) = title(sprintf('%s versus rotation',distance));
      axis([wc(1),wc(end),-inf,inf]); hold on;
      axis('auto y')
      set(th,'fontsize',30);
      FAIRpause;

    else
      figure(k);
      set(vh,'cdata',reshape(Tc,m)')
      setRes(rh,rc);
      subplot(2,2,4);
      set(ph,'visible','off');
      plot(wc(1:j),Dc(1:j),'k-','linewidth',2);
      ph = plot(wc(j),Dc(j),'r.','markersize',20);
      drawnow, 
      FAIRpause(1/100)
    end;
    fprintf('.');
    if ~rem(j,50) || j == length(wc), fprintf('\n'); end;
  end;
  FAIRpause;
end;
%==============================================================================
