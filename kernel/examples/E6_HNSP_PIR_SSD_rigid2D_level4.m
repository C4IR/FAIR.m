%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: Parametric Image Registration
% 
%   - data                 HNSP, Omega=(0,2)x(0,1), level=4, m=[256,128]
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             SSD
%   - pre-registration     rigid2D
%   - optimization         Gauss-Newton
% ===============================================================================

clear, close all, help(mfilename);

setup2DHNSPData; 
imgModel('set','imgModel','splineInter'); 
level = 4; omega = ML{level}.omega; m = ML{level}.m; 
[T,R] = imgModel('coefficients',ML{level}.T,ML{level}.R,omega);
xc = getCellCenteredGrid(omega,m); 
Rc = imgModel(R,omega,xc);

distance('set','distance','SSD');       % initialize distance measure

trafo('reset','trafo','rigid2D');
w0 = trafo('w0'); 

FAIRplots('reset','mode','PIR-rigid','omega',omega,'m',m,'fig',1,'plots',1);
FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m)); 

% ----- call Gaus-Newton ------------------------------------
GNoptn = {'maxIter',500,'Plots',@FAIRplots,'solver','backslash'};
fctn  = @(wc) PIRobjFctn(T,Rc,omega,m,0,[],[],xc,wc);
[wc,his] = GaussNewton(fctn,w0,GNoptn{:});

figure(1); clf
viewImage(imgModel(T,omega,xc),omega,m,'axis','off'); hold on
ph = plotGrid(trafo(wc,xc),omega,m,'spacing',1,'linewidth',1,'color','w');

% plot iteration history
his.str{1} = sprintf('iteration history PIR: distance=%s, y=%s',distance,trafo);
[ph,th] = plotIterationHistory(his,'J',1:4,'fig',2);
%==============================================================================
