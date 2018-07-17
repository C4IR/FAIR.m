%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: SSD and Force Fields
%
%==============================================================================

clear, close all, help(mfilename);
sdiag = @(a) spdiags(reshape(a,[],1),0,length(a),length(a));
setup2DHNSPData; 
level = 7; omega = ML{level}.omega; m = ML{level}.m; % load data
imgModel('reset','imgModel','splineInter')
viewImage('set','axis','off');
[T,R] = imgModel('coefficients',ML{level}.T,ML{level}.R,omega);
X = reshape(getCellCenteredGrid(omega,m),[],2);

[Tc,dT] = imgModel(T,omega,X);
[Rc,dR] = imgModel(R,omega,X);

axT   = [0.55,0.85,0.32,0.47];
axR   = [0.6,0.9,0.20,0.35];
axTF  = [0.55,0.70,0.32,0.395];

FAIRfigure(1,'color','w');
viewImage(Tc,omega,m); hold on; axis off;
ph = plot(axT([1,1,2,2,1]),axT([3,4,4,3,3]),'w-','linewidth',4); 

set(ph,'visible','off');  axis(axT);  

clf; 
viewImage(Rc,omega,m); hold on; axis off;
ph = plot(axR([1,1,2,2,1]),axR([3,4,4,3,3]),'w-','linewidth',4);

set(ph,'visible','off');  axis(axR);  

clf;
viewImage(128+0.5*(Tc-Rc),omega,m); hold on; axis off;
% ph = plot(axT([1,1,2,2,1]),axT([3,4,4,3,3]),'w-','linewidth',4);
% set(ph,'visible','off');  axis(axR);  

X = getCellCenteredGrid(omega,m); X = reshape(X,[],2);
n = prod(m);
gradT     = spdiags(dT,[0,n]);      
maxGradT  = norm(gradT,'inf');
ngradT    = 128+128/maxGradT*gradT;
normGradT = sqrt(sum(gradT.^2,2));
J     = find(normGradT>5e1);
F     = sdiag(Tc-Rc)*gradT;      maxF = norm(F,'inf');
nF    = 128+128/maxF*F;
normF = sqrt(sum(F.^2,2));
K     = find(normF>1e2);

figure(1); clf; set(1,'position',FAIRposition(800),'color','w');
clf; viewImage(ngradT(:,1),omega,m); axis(axT); 
clf; viewImage(ngradT(:,2),omega,m); axis(axT); 

clf; viewImage(128+0.5*Tc,omega,m); hold on;
qh = quiver(X(J,1),X(J,2),gradT(J,1),gradT(J,2),2);
set(qh,'linewidth',3,'color','k'); axis(axTF);

clf; viewImage(128+0.5*Tc,omega,m); hold on;
qh = quiver(X(K,1),X(K,2),F(K,1),F(K,2),2);
set(qh,'linewidth',3,'color','k'); axis(axTF)

%==============================================================================
