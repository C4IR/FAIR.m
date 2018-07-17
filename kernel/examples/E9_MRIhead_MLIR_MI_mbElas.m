%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
% ##2
%==============================================================================
% 
%   - data                 MRI (head), level=7, m=[128,128]
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             MI
%   - pre-registration     affine2D
%   - regularizer          mbElastic
%   - optimization         lBFGS
% ===============================================================================

close all, help(mfilename);

setup2DMRIData
imgModel('reset','imgModel','splineInter','regularizer','moments','theta',1e-2);
distance('reset','distance','MI','nT',8,'nR',8);
trafo('reset','trafo','affine2D');
regularizer('reset','regularizer','mbElastic','alpha',1e-3,'mu',1,'lambda',0);


level = 7; omega = ML{level}.omega; m = ML{level}.m; 

% initialize the interpolation scheme and coefficients
imgModel('reset','imgModel','splineInter'); 
[T,R] = imgModel('coefficients',ML{level}.T,ML{level}.R,omega,'out',0);
xc    = getCellCenteredGrid(omega,m); 
Rc    = imgModel(R,omega,xc);

% initialize distance measure
distance('set','distance','MI','nT',8,'nR',8);       

% initialize regularization, note: yc-yRef is regularized, elastic is staggered 
regularizer('reset','regularizer','mbElastic','alpha',1e-2,'mu',1,'lambda',0);
y0   = getStaggeredGrid(omega,m); yRef = y0; yStop = y0;


% setup and initialize plots 
FAIRplots('set','mode','lBFGS','omega',omega,'m',m,'fig',1,'plots',1);
FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m)); 


% build objective function, note: T coefficients of template, Rc sampled reference
fctn = @(yc) NPIRBFGSobjFctn(T,Rc,omega,m,yRef,yc); fctn([]); % report status


% -- solve the optimization problem -------------------------------------------
[yc,his] = lBFGS(fctn,y0,'maxIter',500,'Plots',@FAIRplots,'yStop',yStop,'tolJ',1e-4);
iter = size(his.his,1)-2; reduction = 100*fctn(yc)/fctn(getStaggeredGrid(omega,m)); yOpt = yc;
fprintf('reduction = %s%% after %d iterations\n',num2str(reduction),iter);
[yc,wc,his] = MLIR(ML,'minLevel',4,'maxLevel',6,'parametric',1,'plotMLiter',0,'plotIter',1);

%==============================================================================
