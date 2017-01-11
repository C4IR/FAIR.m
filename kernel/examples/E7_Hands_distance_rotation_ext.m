%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: various distances versus rotation angle
%
%   - data                 xrays of hands, Omega=(0,20)x(0,25), level=6, m=[128,128]
%   - viewer               viewImage2D
%   - interpolation        splineInter
%   - distance             {'SSD','NCC','MI','NGF'}
%   - transformation       rotation2D
% see also E7_Hands_distance_rotation
%==============================================================================

clear, close all, help(mfilename);

% setup data
setup2DPETCTData; level = 6; omega = ML{level}.omega; m = ML{level}.m; 

DM = {'SSD','NCC','MImex','NGFdot'};


str = @(w) sprintf('T(y(%s^o))',num2str(w*180/pi));
variable = @(k)['DM',DM{k}];

theta = [0,10];

for q=1:length(theta),
    fprintf('============== %s ====================\n\n',...
      sprintf('theta=%s',num2str(theta(q))));
  
    % initialize interpolation [spline] and troansformation [rotation]
    imgModel('reset','imgModel','splineInter','regularizer','moments','theta',theta(q));
    [T,R] = imgModel('coefficients',ML{level}.T,ML{level}.R,omega);
    X  = getCellCenteredGrid(omega,m);
    Rc = imgModel(R,omega,X);
    center = (omega(2:2:end)-omega(1:2:end))'/2;
    trafo('reset','trafo','rotation2D','c',center);
    trafo('w0');

    % run loo over the following distance measures
    edge = 10;
    filename = fullfile(FAIRpath,'temp',...
      sprintf('%s-%s-theta=%s.mat',mfilename,'rot',num2str(theta(q))));
    if ~exist(filename,'file'),
      w = pi/2*linspace(-1,1,101)';
      save(filename,'w','DM');
    else
      clear DM*
      load(filename)
    end;


    for k=1:length(DM),
      var = whos('-file',filename);
      j = find(strcmp({var(:).name},variable(k))==1);

      if isempty(j),
        fprintf('============== %s ====================\n\n',variable(k))
        disp([variable(k),'=dm;']);
        dm = zeros(size(w));
        for j=1:length(w),
          Y = trafo(w(j),X(:));
          Tc = imgModel(T,omega,Y);
          dm(j) = feval(DM{k},Tc,Rc,omega,m,'edge',edge);
          if j == 1,
            figure(k); clf;
            subplot(1,3,1); viewImage(Rc,omega,m); title('reference');
            subplot(1,3,2); ph = viewImage(Tc,omega,m); th = title(str(w(j)));
            subplot(1,3,3); rh = plot(w(j),dm(j),'r.','markersize',20);
            axis([w(1),w(end),-inf,inf]); hold on;      title(DM{k})
            axis('auto y')
          else
            set(ph,'cdata',reshape(Tc,m)'); set(th,'string',str(w(j)));
            subplot(1,3,3); set(rh,'visible','off');
            plot(w(1:j),dm(1:j),'k-','linewidth',2);
            rh = plot(w(j),dm(j),'r.','markersize',20);    pause(1/100)
          end;
          fprintf('.'); if ~rem(j,50) || j==length(w), fprintf('\n'); end;
        end;
        eval([variable(k),'=dm;']);
        save(filename,'-append',variable(k));
        fprintf('============== %s done ===============\n\n',variable(k))
      end;
    end;

    Name = @(k) sprintf('%s-theta=%s',DM{k},num2str(theta(q)));
    for k=1:4,
      eval(['dm=',variable(k),';']);
      dm = dm(1:length(w));
      [wOpt,j] = min(dm);
      FAIRfigure(k);
      ph = plot(w,dm,'-',w(j),dm(j),'*');
      set(ph,'linewidth',2,'color','k','markersize',20);
      set(gca,'fontsize',30);
      a = max(dm)-0.2*(max(dm)-min(dm));
      th = text(w(j),a,['$w^*=',num2str(w(j)),'$']);
      set(th,'fontsize',30,'interpreter','latex','horizontalalignment','center');
    end;
    fprintf('=====================================================================\n');
    
end;
return
% =========================================================================
