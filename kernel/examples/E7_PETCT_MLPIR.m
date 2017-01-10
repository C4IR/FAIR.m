%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: various distances and Multi-Level Parametric Image Registration
%
%   - data                 PETCT, Omega=(0,140)x(0,151), level=4:7, m=[128,128]
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             {'SSD','NCC','MI','NGF'}
%   - transformation       affine2D
% see also E7_PETCT_MLPIR_ext
%==============================================================================

clear, close all, help(mfilename);
setup2DPETCTData;                                         % load data

% a list of distance measures to be used
DM = {'SSD','NCC','MI','NGF'};

OPTpara = FAIRcell2struct(optPara('PIR-GN'));

for dm = 1:length(DM), % run over all distance measures

  % initialize interpolation, using a smooth representation (theta=1e0) 
  imgModel('reset','imgModel','splineInter','regularizer','moments','theta',1e0);

  % initialize transformation, create initial guess and reference for stopping
  trafo('reset','trafo','affine2D'); wStop = trafo('w0'); w0 = wStop;

  % initialize distance and display options
  distance('reset','distance',DM{dm}); distance('disp')
  
  % run MLPIR using sufficient amount of details (level=5)
  wSmooth =  MLPIR(ML,'minLevel',5,'plotIter',0,'plotMLiter',0);

  % refine interpolation (theta=1e-3)
  imgModel('reset','imgModel','splineInter','regularizer','moments','theta',1e-3);
  level = length(ML); omega = ML{level}.omega; m = ML{level}.m;
  [T,R] = imgModel('coefficients',ML{level}.T,ML{level}.R,omega);
  
  % start PIR, using the result from the smooth problem as starting guess
  
  % initialize plots
  FAIRplots('set','mode','PIR');
  FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m));

  % optimize
  xc = getCellCenteredGrid(omega,m);   
  Rc = imgModel(R,omega,xc);
  fctn = @(wc) PIRobjFctn(T,Rc,omega,m,0,[],[],xc,wc); fctn([]);
  [wc,his] = GaussNewton(fctn,w0,OPTpara{:});

  % visualize results
  yc = trafo(wc,xc);
  R0 = imgModel(R,omega,xc);
  T0 = imgModel(T,omega,xc);
  Tc = imgModel(T,omega,yc);

  figure(11); clf;
  viewImage(T0,omega,m,'axis','off'); hold on;
  plotGrid(yc,omega,m,'spacing',ceil(m/32),'linewidth',2,'color','w');

  figure(12); clf;
  overlayImage2D(Tc,R0,omega,m); axis off;
end;

%==============================================================================
