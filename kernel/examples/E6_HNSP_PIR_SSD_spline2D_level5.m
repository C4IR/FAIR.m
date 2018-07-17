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
%   - data                 HNSP, Omega=(0,2)x(0,1), level=5, m=[256,128]
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             SSD
%   - pre-registration     splineTransformation2D
%   - optimization         Gauss-Newton
% 
% ===============================================================================

clear, close all, help(mfilename);

setup2DHNSPData; 
imgModel('set','imgModel','splineInter'); 
level = 5; omega = ML{level}.omega; m = ML{level}.m; 
[T,R] = imgModel('coefficients',ML{level}.T,ML{level}.R,omega);
xc = getCellCenteredGrid(omega,m); 
Rc = imgModel(R,omega,xc);

distance('set','distance','SSD');       % initialize distance measure

% initialize transformation and starting guess; 
% here: spline with 2 times [4,5] coefficients
trafo('reset','trafo','splineTransformation2D','omega',omega,'m',m,'p',[4,5]);
w0 = trafo('w0'); 

FAIRplots('reset','mode','PIR-spline','omega',omega,'m',m,'fig',1,'plots',1);
FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m)); 

% ----- call Gaus-Newton ------------------------------------
GNoptn = {'maxIter',50,'Plots',@FAIRplots,'solver','backslash'};
fctn  = @(wc) PIRobjFctn(T,Rc,omega,m,0,[],[],xc,wc);
[wc,his] = GaussNewton(fctn,w0,GNoptn{:});

figure(1); clf
viewImage(imgModel(T,omega,xc),omega,m,'axis','off'); hold on
ph = plotGrid(trafo(wc,xc),omega,m,'spacing',1,'linewidth',1,'color','w');
%==============================================================================
