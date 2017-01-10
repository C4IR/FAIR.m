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
% see also E7_PETCT_MLPIR
%==============================================================================

clear, close all, help(mfilename);

setup2DPETCTData;

FAIRfigure(12,'position','default');
FAIRfigure(13,'position','default');
FAIRfigure(14,'position','default');

FAIRprint = @(varargin) [];
FAIRprint('reset','folder',fullfile(FAIRpath,'../temp','PETCT'),...
  'pause',1,'obj','gca','format','jpg','draft','off');
FAIRprint('disp');

rescaleGrayvalues = @(T0) (255*(T0-min(T0))/max(T0-min(T0))).^3;

mname = mfilename;
File  = @(str) sprintf('%s-%s',mname,str);

outpath = fullfile(FAIRpath,'temp','PETCT');
if ~exist(outpath,'dir'),
  mkdir(outpath);
end;

Write = @(T,str) imwrite(uint8(flipud(reshape(T,m)')),...
  fullfile(outpath,[File(str),'.jpg']));

DM = {'SSD','NCC','MImex','NGFdot'};
map = @(T) 160*(T/160).^(0.5);

for dm = 1:length(DM),

  imgModel('reset','imgModel','splineInter','regularizer','moments','theta',1e0);
  level = 5; omega = ML{level}.omega; m = ML{level}.m;
  [T,R] = imgModel('coefficients',ML{level}.T,ML{level}.R,omega);
  trafo('reset','trafo','affine2D');  % initialize transformation
  wStop = trafo('w0');
  wOpt = wStop;

  distance('reset','distance',DM{dm});
  distance('disp');
  
  if 1,  
    wSmooth =  MLPIR(ML,'minLevel',5,'plotIter',0,'plotMLiter',0);

    imgModel('reset','imgModel','splineInter','regularizer','moments','theta',1e-3);
    level = length(ML); omega = ML{level}.omega; m = ML{level}.m;
    [T,R] = imgModel('coefficients',ML{level}.T,ML{level}.R,omega);
    xc = getCellCenteredGrid(omega,m);
    Rc = imgModel(R,omega,xc);
    fctn = @(wc) PIRobjFctn(T,Rc,omega,m,0,[],[],xc,wc);
    wOpt = GaussNewton(fctn,wSmooth,'yStop',wStop);
  end;

  omega = ML{end}.omega; m = ML{end}.m;
  xc = getCellCenteredGrid(omega,m);
  R0 = imgModel(R,omega,xc);
  T0 = imgModel(T,omega,xc);

  if dm == 0,
    FAIRfigure(13)
    str = 'R0'; var = R0;
    viewImage(var,omega,m,'axis','off','title',File(str));  Write(var,str);

    FAIRfigure(12)
    str = 'T0'; var = T0;
    viewImage(var,omega,m,'axis','off','title',File(str));  Write(var,str);

    FAIRfigure(14); clf;
    Tk = rescaleGrayvalues(T0);
    overlayImage2D(Tk,R0,omega,m);
    FAIRprint(File(['D0']));
    return
  end;

  Yopt = trafo(wOpt,xc);
  Tc = imgModel(T,omega,Yopt);

  FAIRfigure(12); clf;
  str = ['G-',distance];  var = T0;
  viewImage(var,omega,m,'axis','off'); hold on;
  plotGrid(grid2grid(Yopt,m,'centered','nodal'),omega,m,...
    'spacing',ceil(m/16),'linewidth',3,'color','w');
  set(gca,'position',[0 0 1 1]);
  FAIRprint(File(['G-',distance]));

  FAIRfigure(13); clf;
  str = ['T-',distance]; var = Tc;
  viewImage(var,omega,m,'axis','off','title',File(str));  Write(var,str);

  FAIRfigure(14); clf;
  Tk = rescaleGrayvalues(Tc);
  overlayImage2D(Tk,R0,omega,m);
  FAIRprint(File(['D-',distance]));
end;

FAIRpause(5); close all

%==============================================================================
