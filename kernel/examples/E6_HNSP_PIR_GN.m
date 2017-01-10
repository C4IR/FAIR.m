%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: MLPIR, Parametric Image Registration
% 
%   - data                 HNSP, Omega=(0,2)x(0,1), level=4, m=[256,128]
%   - viewer               viewImage2D
%   - interpolation        splineInterMex
%   - distance             SSD
%   - pre-registration     affine2D
%   - optimization         Gauss-Newton
% ===============================================================================

clear; help(mfilename);

% setup data (MultiLevel based) and initialize image viewer
setup2DHNSPData; 

% initialize the interpolation scheme and its coefficients
imgModel('set','imgModel','splineInterMex'); 
level = 8; omega = ML{level}.omega; m = ML{level}.m; 
[T,R] = imgModel('coefficients',ML{level}.T,ML{level}.R,omega);
xc    = getCellCenteredGrid(omega,m); 
Rc    = imgModel(R,omega,xc);

% initialize distance measure
distance('set','distance','SSD');       

% initialize the transformation and a starting guess
trafo('reset','trafo','affine2D');
w0 = trafo('w0'); 

% setup plots and initialize
FAIRplots('reset','mode','PIR-Gauss-Newton','omega',omega,'m',m,'fig',1,'plots',1);
FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m)); 

% build objective function
% note: T  is template image
%       Rc is sampled reference
%       optional Tikhonov-regularization is disabled by setting m = [], wRef = []
%       beta = 0 disables regularization of Hessian approximation
beta = 0; M = []; wRef = [];
fctn = @(wc) PIRobjFctn(T,Rc,omega,m,beta,M,wRef,xc,wc); 
fctn([]);   % report status

% -- solve the optimization problem -------------------------------------------
[wc,his] = GaussNewton(fctn,w0,'Plots',@FAIRplots,'solver',[],'maxIter',100);
plotIterationHistory(his,'J',[1,2,5],'fig',20+level); 
%==============================================================================
