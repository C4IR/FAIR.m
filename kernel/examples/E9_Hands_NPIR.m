%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
% ##2
%==============================================================================
%
% Example for single level image registration of hand data
%
%   - data                 Hands, Omega=(0,20)x(0,25), level=5, m=[32,32]
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             SSD
%   - regularizer          mbElastic
%   - optimization         Gauss-Newton
% ===============================================================================

% setup data, pick level 5, and initialize image viewer
setup2DhandData;
level = 4;
omega = ML{level}.omega; 
m     = ML{level}.m; 

% initialize the image model, distance measure and regularizer 
imgModel('reset','imgModel','splineInter'); 
[T,R] = imgModel('coefficients',ML{level}.T,ML{level}.R,omega,'out',0);
xc    = getCellCenteredGrid(omega,m); 
Rc    = imgModel(R,omega,xc);
distance('set','distance','SSD');       
regularizer('reset','regularizer','mfElastic','alpha',1e4,'mu',1,'lambda',0);
y0   = getStaggeredGrid(omega,m); yRef = y0; yStop = y0;

% build objective function, note: T coefficients of template, Rc sampled reference
fctn = @(yc) NPIRobjFctn(T,Rc,omega,m,yRef,yc); 
fctn([]); % report status of opjective function

% %% run default: Gauss-Newton with multigrid solver
NPIRpara = FAIRcell2struct(optPara('NPIR-GN'));

FAIRplots('reset','mode','Hands-GaussNewton-MG','fig',1);
FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m)); 
yOpt =  GaussNewton(fctn,y0,NPIRpara{:});


%==============================================================================
