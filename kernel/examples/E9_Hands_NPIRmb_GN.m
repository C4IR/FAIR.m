%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
% ##2
%==============================================================================
% 
%   - data                 Hand, Omega=(0,20)x(0,25), 
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             MI
%   - regularizer          mbElastic
%   - optimizer            Gauss-Newton
% ===============================================================================

clear; close all; help(mfilename)

setup2DhandData
level = 5; omega = ML{level}.omega; m = ML{level}.m;
imgModel('reset','imgModel','splineInter','regularizer','moments','theta',0.01);
distance('reset','distance','MI','nT',10,'nR',10);
regularizer('reset','regularizer','mbElastic','alpha',1000,'mu',1,'lambda',0);

[T,R] = imgModel('coefficients',ML{level}.T,ML{level}.R,omega);
xc    = getStaggeredGrid(omega,m); % starting guess and reference for regularization
Rc    = imgModel(R,omega,center(xc,m)); 

% - initialize FAIR plots
FAIRplots('set','mode','NPIR-mb','fig',1);
FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m)); 

% - run Non-Parametric Image Registration (Gauss-Newton)
fctn = @(yc) NPIRobjFctn(T,Rc,omega,m,xc,yc);  fctn([]);
NPIRpara = FAIRcell2struct(optPara('NPIR-GN','Plots',@FAIRplots));
yc = GaussNewton(fctn,xc,NPIRpara{:});
%==============================================================================
