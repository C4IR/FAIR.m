%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: Nesterov,  Parametric Image Registration
% 
%   - data                 Hand, Omega=(0,20)x(0,25), level=5, m=[32,32]
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             SSD
%   - transformation       affine2D
%   - optimization         Nesterov
% ===============================================================================

clear, close all, help(mfilename);

setup2DhandData;
viewImage('reset','viewImage','viewImage2D','colormap',bone(256),'axis','off');
imgModel('reset','imgModel','splineInter','regularizer','moments','theta',1e-1);

level = 5; 
omega = ML{level}.omega; 
m     = ML{level}.m;

[T,R] = imgModel('coefficients',ML{level}.T,ML{level}.R,omega,'out',0);

distance('reset','distance','SSD');

trafo('reset','trafo','affine2D'); 
w0 = trafo('w0'); beta = 0; M =[]; wRef = []; % disable regularization

% initialize plots
FAIRplots('reset','mode','PIR-Nesterov','fig',1);
FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m)); 

xc = getCellCenteredGrid(omega,m); 
Rc = imgModel(R,omega,xc);
fctn = @(wc) PIRobjFctn(T,Rc,omega,m,beta,M,wRef,xc,wc);

[Jc,para,dJ] = fctn(w0); 
stepLength = 0.5*norm(dJ);

wc = Nesterov(fctn,w0,'maxIter',150,'lipschitzC',stepLength);

%==============================================================================
