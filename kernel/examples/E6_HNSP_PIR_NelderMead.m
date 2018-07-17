%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: PIR,  Parametric Image Registration
% 
%   - data                 Hand, Omega=(0,20)x(0,25), level=5, m=[32,32]
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             SSD
%   - pre-registration     affine2D
%   - optimization         Nelder-Mead (fminsearch)
% ===============================================================================

clear, close all, help(mfilename);

% setup data, interpolation, transformation, distance, 
% extract data from a particular level amd compute coefficients for interpolation
setup2DHNSPData; 
imgModel('set','imgModel','splineInter'); 
level = 4; omega = ML{level}.omega; m = ML{level}.m; 
[T,R] = imgModel('coefficients',ML{level}.T,ML{level}.R,omega);
xc = getCellCenteredGrid(omega,m); 
Rc = imgModel(R,omega,xc);
distance('set','distance','SSD');       
trafo('reset','trafo','affine2D');
w0 = trafo('w0'); 

% build objective function
% note: T  is template image
%       Rc is sampled reference
%       optional Tikhonov-regularization is disabled by setting m = [], wRef = []
%       beta = 0 disables regularization of Hessian approximation
beta = 0; M = []; wRef = [];
fctn = @(wc) PIRobjFctn(T,Rc,omega,m,beta,M,wRef,xc,wc); 
fctn([]);   % report status

optn = optimset('display','iter');
[ycNM,JcNM,exitFlag,out] = fminsearch(fctn,w0,optn);
yc = trafo(ycNM,xc);
showResults(ML,yc,'level',level)
%==============================================================================
