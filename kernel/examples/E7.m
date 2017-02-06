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

function E7(example)

if nargin == 0, example = 'PETCT';  end;
imagePath  = fullfile(FAIRpath,'temp',example);
switch example
  case 'HNSP'
    setup2DHNSPData;
    theta = 1;
  case 'MRIhead',
    setup3DMRIData;
    theta = 1;
  case 'PETCT',
    setup2DPETCTData;
    theta = 0;
  otherwise, return;
end
level = 6; omega = ML{level}.omega; m = ML{level}.m;

DM = {'SSD','NCC','MImex','NGFdot'};
viewImage('set',viewPara{:})
filename = fullfile(FAIRpath,'temp',sprintf('%s-%s-%s.mat',mfilename,'rotation',example));

% 2e1 MRIhead
imgModel('reset','imgModel','splineInter','regularizer','moments','theta',theta);
[T,R] = imgModel('coefficients',ML{level}.T,ML{level}.R,omega);
xc    = getCellCenteredGrid(omega,m); 
xc    = reshape(xc,[],2);
Tin   = imgModel(T,omega,xc);
Rin   = imgModel(R,omega,xc);

%%  plots
k = figureh(1); clf; set(k,'color','w');
subplot(1,2,1); viewImage(Tin,omega,m);
subplot(1,2,2); viewImage(Rin,omega,m);
FAIRpause(1);

center = (omega(2:2:end)-omega(1:2:end))'/2;
trafo('reset','trafo','rotation2D','c',center);
trafo('w0');
if ~exist(filename,'file'), 
  w = pi/4*linspace(-1,1,101)';
  save(filename,'w','DM');  
else
  clear DM*
  load(filename)
end;

edge = 100;

for k=1:length(DM),
  variable = ['DM',DM{k}];
  var = whos('-file',filename);
  j = find(strcmp({var(:).name},variable)==1);

  if isempty(j),
    fprintf('============== %s ====================\n\n',variable)
    disp([variable,'=dm;']);
    dm = zeros(size(w));
    for j=1:length(w),
      Y  = trafo(w(j),xc(:));
      Tc = imgModel(T,omega,Y);
      dm(j) = feval(DM{k},Tc,Rin,omega,m,'edge',edge);
      if j== 1,
        figureh(3);
        ph = viewImage(Tc,omega,m);
      else
        set(ph,'cdata',reshape(Tc,m)'); drawnow
        FAIRpause(1/100)
      end;
      title(sprintf('%d/%d:%d/%d',k,length(DM),j,numel(w)));
    end;
    eval([variable,'=dm;']);
    save(filename,'-append',variable);
  end;
end;

%%
load(filename)


filename = fullfile(FAIRpath,'temp',...
  sprintf('%s-%s-%s.mat',mfilename,'rotation',example));



for k=1:4,
  variable = ['DM',DM{k}];
  eval(['dm=',variable,';']);
  dm = dm(1:length(w));
  [wOpt,j] = min(dm);
  figureh(k); clf; set(k,'color','w');
  ph = plot(w,dm,'-',w(j),dm(j),'*');
  set(ph,'linewidth',2,'color','k','markersize',20);
  set(gca,'fontsize',30);
  a = max(dm)-0.2*(max(dm)-min(dm));
  th = text(w(j),a,['$w^*=',num2str(w(j)),'$']);
  set(th,'fontsize',30,'interpreter','latex','horizontalalignment','center');
end;


%% finish
filename = fullfile(FAIRpath,'temp',...
  sprintf('%s-%s-%s.mat',mfilename,'translation',example));


trafo('reset','trafo','translation2D');
if ~exist(filename,'file'), 
  [w1,w2] = ndgrid(0.1*(omega(2)-omega(1))*linspace(-1,1,21),...
                   0.2*(omega(4)-omega(3))*linspace(-1,1,21));
  save(filename,'w1','w2','DM');  
else
  clear DM*
  load(filename)
end;

for k=1:length(DM),
  variable = ['DM',DM{k}];
  disp([variable,'=dm;']);

  var = whos('-file',filename);
  j = find(strcmp({var(:).name},variable)==1);

  if isempty(j),
    fprintf('============== %s ====================\n\n',variable)
    dm = zeros(size(w1));
    for j=1:numel(w1),
      Y = trafo([w1(j);w2(j)],xc(:));
      Tc = imgModel(T,omega,Y);
      dm(j) = feval(DM{k},Tc,Rin,omega,m,'edge',edge);
      if j== 1,
        figureh(3);
        ph = viewImage(Tc,omega,m);
      else
        set(ph,'cdata',reshape(Tc,m)'); drawnow
        pause(1/100)
      end;
      title(sprintf('%d/%d:%d/%d',k,length(DM),j,numel(w1)));
    end;
    eval([variable,'=dm;']);
    save(filename,'-append',variable);
  end;
end;

load(filename)

%figureh(1); close(1); figureh(1); clf; set(1,'position',position(800),'color','w');
%figureh(2); close(2); figureh(2); clf; set(2,'position',position(800),'color','w');
shift = [3000,0.5,-1,0.08];
Ztick = {[0:1000:4000],[0:0.1:0.3],[-1:0:1],[0:0.02:0.06]};

for k=1:4,

  variable = ['DM',DM{k}];
  eval(['dm=',variable,';']);

  fig = figureh(10+k);
  if ~isnumeric(fig), fig = fig.Number; end;
  close(fig); figureh(fig); set(fig,'color','w');
  ph=mesh(w1,w2,dm-shift(k)); grid off;
  set(gca,'fontsize',30);
  view(-135,25);
  set(ph,'linewidth',2)
  
  fig = figureh(20+k);
  if ~isnumeric(fig), fig = fig.Number; end;
  close(fig); figureh(fig); clf; set(fig,'color','w');
  contour(w1,w2,dm,10,'linewidth',2);
  set(gca,'fontsize',30);
end;
%==============================================================================
