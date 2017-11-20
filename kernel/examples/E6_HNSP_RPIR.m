%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: Regularized Parametric Image Registration,  by spline-moment matrix
% 
%   - data                 HNSP, Omega=(0,2)x(0,1), level=5, m=[256,128]
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             SSD
%   - pre-registration     splineTransformation2D, regularized by spline-moment matrix
%   - optimization         Gauss-Newton
% ===============================================================================

clear, close all, help(mfilename);
%%
setup2DHNSPData;  level = 5;  omega = ML{level}.omega; m = ML{level}.m; 
p = [8,8];
imgModel('reset','imgModel','splineInter','regularizer','none','theta',0);
[T,R] = imgModel('coefficients',ML{level}.T,ML{level}.R,omega);
distance('reset','distance','SSD');
trafo('reset','trafo','splineTransformation2D','omega',omega,'m',m,'p',p);
w0 = trafo('w0');                        % get starting guess and stopping

% update the regularizer for parametric image registration
hd = prod((omega(2:2:end)-omega(1:2:end))./m);
Mi = @(i) 1e2*toeplitz([96,-54,0,6,zeros(1,p(i)-4)]);
Qi = @(i) toeplitz([120.8,59.55,6,0.05,zeros(1,p(i)-4)])/7;
M1 = kron(speye(p(2)),sparse(Mi(1)));
M2 = kron(sparse(Mi(2)),speye(p(1)));
M  = hd*sparse(kron(speye(2),...
  kron(Qi(2),Mi(1))+2*kron(Mi(2),Mi(1))+kron(Mi(2),Qi(1))));

%% set up elastic regularization matrix
mu     = 2;
lambda = 1;
B = getElasticMatrixNodal(omega,m,mu,lambda);
A = B'*B;
[y,dy] = splineTransformation2D(w0,getNodalGrid(omega,m),'m',m+1);
M = [dy'*A*dy];

%%
% M = speye(length(w0));
alpha = 10;
beta  = 0;
wRef  = w0;

% initialize objective function
xc = getCellCenteredGrid(omega,m);
Rc = imgModel(R,omega,xc);
fctn  = @(wc) PIRobjFctn(T,Rc,omega,m,beta,alpha*M,wRef,xc,wc);
fctn([]);

% initialize plots
FAIRplots('reset','mode','PIR-regularized','omega',omega,'m',m,'fig',1,'plots',1);
FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m)); 

xc = getCellCenteredGrid(omega,m); 
Rc = imgModel(R,omega,xc);

% ----- call Gaus-Newton ------------------------------------
GNoptn = {'maxIter',100,'tolY',1e-4,'tolJ',1e-4,'tolG',1,'solver','backslash'};
[wOpt,his] = GaussNewton(fctn,w0,GNoptn{:},'Plots',@FAIRplots);

% plot iteration history
his.str{1} = sprintf('iteration history PIR: distance=%s, y=%s',distance,trafo);
[ph,th] = plotIterationHistory(his,'J',1:4);
%==============================================================================
