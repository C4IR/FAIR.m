%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% 
%   - data                 3d-brains, Omega=(0,20)x(0,10)x(0,20), level=3:6, m=[128,64,128]
%   - viewer               imgmontage
%   - interpolation        linearInterMex
%   - distance             SSD
%   - pre-registration     affine3D
%   - regularizer          mfElastic
%   - optimizer            Gauss-Newton
% ===============================================================================

setup3DbrainData;

% extract data of level 4
level = 4; omega = ML{level}.omega; m = ML{level}.m; 
imgModel('reset','imgModel','linearInterMex'); 
[T,R] = imgModel('coefficients',ML{level}.T,ML{level}.R,omega);
xc = getCellCenteredGrid(omega,m); 
Rc = imgModel(R,omega,xc);

% initialize distance measure
distance('reset','distance','SSD');       

% initialize the transformation and a starting guess
% trafo('reset','trafo','affine3Dsparse');
trafo('reset','trafo','splineTransformation3Dsparse',...
  'omega',omega,'p',[ 3 4 5],'m',m);
w0 = trafo('w0');


% setup plots and initialize
FAIRplots('reset','mode','PIR-GN','omega',omega,'m',m,'fig',1,'plots',1);
FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m)); 

%------------------------------------------------------------------------------
% build objective function
% note: T  is data for template image
%       Rc is sampled reference image
%       optional Tikhonov-regularization is disabled by setting m = [], wRef = []
%       beta = 0, M = [], wRef = []:  
%       disables additional regularization of Hessian approximation
beta = 0; M = []; wRef = [];
fctn = @(wc) PIRobjFctn(T,Rc,omega,m,beta,M,wRef,xc,wc); 
fctn([]);   % report status

%------------------------------------------------------------------------------
%% -- solve the optimization problem on one level
[wc,his] = GaussNewton(fctn,w0,'solver','CG','Plots',@FAIRplots); return;

%------------------------------------------------------------------------------
%% finally: run the MultiLevel Non-Parametric Image Registration
[wc,his] = MLPIR(ML);

%==============================================================================


