%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: distances and SSD forces for hand data
%
%   - data                 Hand, Omega=(0,20)x(0,25), level=3:7, m=[128,128]
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             MI
%   - transformation       translation2D
%==============================================================================

clear, close all, help(mfilename);

fprintf('setup data, viewer, interpolation, transformation\n');
setup2DhandData; level = 6; omega = ML{level}.omega; m = ML{level}.m;
viewImage('reset',viewPara{:},'axis','off');
imgModel('reset','imgModel','splineInter');
trafo('reset','trafo','translation2D');
R  = imgModel('coefficients',ML{level}.R,[],omega,'out',0);
xc = reshape(getCellCenteredGrid(omega,m),[],2);
yc = reshape(trafo([0.75;0],xc(:)),[],2);

fprintf('%20s : %s\n','viewImage',viewImage);
fprintf('%20s : %s\n','imgModel',imgModel);

% compute R(xc), T(yc), dT and forces
Rc = imgModel(R,omega,xc);
[Tc,dT] = imgModel(R,omega,yc);
dT = full(spdiags(dT,[0,size(dT,1)]));
F  = spdiags(Tc-Rc,0,length(Tc),length(Tc))*dT;
  
% compute lengthes for plots
fac = 128/norm(dT(:),'inf');
lengthdT = sqrt(sum(dT.^2,2));
normdT   = max(lengthdT);
J  = find(lengthdT>1e-2*normdT);
lengthF = sum(F.^2,2);
normF   = sqrt(max(lengthF));
K  = find(lengthF>1e-2*normF);
  
% plot R, T, T-R
FAIRfigure(1,'figname',mfilename); clf;
subplot(1,3,1);  viewImage(Rc,omega,m);            th(1) = title('R');
subplot(1,3,2);  viewImage(Tc,omega,m);            th(2) = title('T');
subplot(1,3,3);  viewImage(128+(Tc-Rc)/2,omega,m); th(3) = title('T-R');
set(th,'fontsize',30);
FAIRpause;
  
% plot \partial_j T and \nabla T
figure(2); clf;  colordef(gcf,'black');
subplot(1,3,1);  viewImage(128+fac*dT(:,1),omega,m);    th(1) = title('\partial_1T');
subplot(1,3,2);  viewImage(128+fac*dT(:,2),omega,m);    th(2) = title('\partial_2T');
subplot(1,3,3);  viewImage(128+Tc/2,omega,m); hold on;
qh = quiver(xc(J,1),xc(J,2),dT(J,1)/normdT,dT(J,2)/normdT,2);
set(qh,'color','b','linewidth',1.5);
th(3) = title('\nabla T');
set(th,'fontsize',30);
FAIRpause;
  
% plot \nabla T, T-R, forces
figure(3); clf;  colordef(gcf,'black');
subplot(1,3,1);  viewImage(128+Tc/2,omega,m); hold on;
qh = quiver(xc(J,1),xc(J,2),dT(J,1)/normdT,dT(J,2)/normdT,0);
set(qh,'color','b','linewidth',1.5);
th(1) = title('\nabla T');

subplot(1,3,2);  vh = viewImage(128+(Tc-Rc)/2,omega,m); th(2) = title('T-R');
  
subplot(1,3,3);  viewImage(128+Tc/2,omega,m); hold on;
qh = quiver(xc(K,1),xc(K,2),F(K,1)/normF,F(K,2)/normF,2);
set(qh,'color','r','linewidth',1.5);
th(3) = title('forces');
set(th','fontsize',30);

%==============================================================================
