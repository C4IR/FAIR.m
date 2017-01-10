%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
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
%   - pre-registration     affine2D
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
fctn = @(yc) NPIRobjFctn(T,Rc,omega,m,yRef,yc); fctn([]); % report status


% %% run default: Gauss-Newton with multigrid solver
% NPIRopt = optPara('NPIR-GN');
% NO = FAIRcell2struct(NPIRopt);
% FAIRplots('reset','mode','Hands-GaussNewton-MG','fig',1);
% FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m)); 
% yOpt =  GaussNewton(fctn,y0,NO{:});
% 
% FAIRpause(2);

%% run over a variety of optimizer

[Jc,para,dJ] = fctn(y0); 
lipschitzC   = 0.5*norm(dJ);
stepLength   = 0.5/norm(dJ);

% first order, no solvers
methods(1).name    = 'Hands-SteepestDescent';
methods(1).scheme  = @SteepestDescent;
methods(1).para    = FAIRcell2struct(optPara('SD',...
  'stepLength',stepLength,'maxIter',50));

methods(2).name    = 'Hands-Nesterov';
methods(2).scheme  = @Nesterov;
methods(2).para    = FAIRcell2struct(optPara('Nesterov','lipschitzC',lipschitzC));

% Nesterov with Gauss-Newtom
methods(3).name    = 'Hands-Nesterov-GN';
methods(3).scheme  = @NesterovGN;
methods(3).para    = FAIRcell2struct(optPara('Nesterov',...
  'lipschitzC',lipschitzC,'solver','MG-elastic','maxIter',10));

% Gauss-Newton style
methods(4).name    = 'Hands-GaussNewton-MG';
methods(4).scheme  = @GaussNewton;
methods(4).para    = FAIRcell2struct(optPara('GN','solver','MG-elastic'));

methods(5).name    = 'Hands-GaussNewton-mf-CG-elastic';
methods(5).scheme  = @GaussNewton;
methods(5).para    = FAIRcell2struct(optPara('GN','solver','CG-elastic'));

methods(6).name    = 'Hands-GaussNewton-mf-PCG-elastic';
methods(6).scheme  = @GaussNewton;
methods(6).para    = FAIRcell2struct(optPara('GN','solver','PCG-elastic'));

methods(7).name    = 'Hands-GaussNewton-backslash';
methods(7).scheme  = @GaussNewton;
methods(7).para    = FAIRcell2struct(optPara('GN','solver','backslash'));

methods(8).name    = 'Hands-GaussNewton-PCG-SGS';
methods(8).scheme  = @GaussNewton;
methods(8).para    = FAIRcell2struct(optPara('GN','solver','mbPCG-SGS'));

methods(9).name    = 'Hands-GaussNewton-PCG-ICHOL';
methods(9).scheme  = @GaussNewton;
methods(9).para    = FAIRcell2struct(optPara('GN','solver','mbPCG-ICHOL'));

methods(10).name   = 'Hands-GaussNewton-PCG-Jacobi';
methods(10).scheme = @GaussNewton;
methods(10).para   = FAIRcell2struct(optPara('GN','solver','mbPCG-Jacobi'));

methods(11).name   = 'Hands-GaussNewton-CG';
methods(11).scheme = @GaussNewton;
methods(11).para   = FAIRcell2struct(optPara('GN','solver','mbCG'));

methods(12).name   = 'Hands-lBFGS';
methods(12).scheme = @lBFGS;
methods(12).para   = FAIRcell2struct(optPara('lBFGS','solver','backslash'));

methods(13).name   = 'Hands-TrustRegion';
methods(13).scheme = @TrustRegion;
methods(13).para   = FAIRcell2struct(optPara('TrustRegion','solver','backslash'));

K = setdiff(1:length(methods),3);
% K(1:4) = []
for k=K
  a = find(strcmp(methods(k).para,'solverr'));
  if isempty(a),
    solver = 'is empty';
  else
    solver = methods(k).para(a+1);
  end;
  
  fprintf('test %s // solver %s\n',methods(k).name,solver);
  
  if ismember(k,[1:6,13]), % matrix free code
    regularizer('reset','regularizer','mfElastic','alpha',1e4,'mu',1,'lambda',0);
  else
    regularizer('reset','regularizer','mbElastic','alpha',1e4,'mu',1,'lambda',0);
  end;
  FAIRmessage(func2str(methods(k).scheme),'#')
  
  % - initialize FAIR plots
  FAIRplots('set','mode',func2str(methods(k).scheme),'fig',1);
  FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m));

  yOpt =  methods(k).scheme(fctn,y0,methods(k).para{:});
  FAIRpause(2);
end
%==============================================================================
