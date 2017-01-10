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
% see also E7_Hands_distance_rotation_ext
%==============================================================================

clear, close all, help(mfilename);

% setup data
setup2DhandData; 
level = 6; omega = ML{level}.omega; m = ML{level}.m; 

DM = {'SSD','NCC','MI','NGF'}; % the distance measures
str = @(w) sprintf('T(y(%s^o))',num2str(w*180/pi)); % used in plots for titles

% initialize interpolation, here spline interpolation with various theta's 
theta= 0;                                           % play with theta!
imgModel('reset','imgModel','splineInter','regularizer','moments','theta',theta);
[T,R] = imgModel('coefficients',ML{level}.T,ML{level}.R,omega);

% initialize a grid xc and compute Rc = R(xc)
xc  = getCellCenteredGrid(omega,m);
Rc = imgModel(R,omega,xc);

% initialize the transformation, here rotation in 2D
center = (omega(2:2:end)-omega(1:2:end))'/2;
trafo('reset','trafo','rotation2D','c',center);

% parameter runs in [-pi/2,pi/2]
wc = pi/2*linspace(-1,1,101)';

% run loop over the following distance measures

for k=1:length(DM),
  fprintf('============== %s ========================\n',DM{k})
  dc = zeros(size(wc));                  % allocate memory for D(w(j))
  for j=1:length(wc),                    % loop over all parameters w(j)    
    yc = trafo(wc(j),xc(:));             % compute the transformation 
    Tc = imgModel(T,omega,yc);              % compute transformed template
    dc(j) = feval(DM{k},Tc,Rc,omega,m);  % evaluate distance measure

    % do some plots
    if j == 1,                          % initialize the plot
      figure(k); clf;
      subplot(1,3,1); viewImage(Rc,omega,m); title('reference');
      subplot(1,3,2); ph = viewImage(Tc,omega,m); th = title(str(wc(j)));
      subplot(1,3,3); rh = plot(wc(j),dc(j),'r.','markersize',20);
      axis([wc(1),wc(end),-inf,inf]); hold on;      title(DM{k})
      axis('auto y')
    else                                % update the plot
      set(ph,'cdata',reshape(Tc,m)'); set(th,'string',str(wc(j)));
      subplot(1,3,3); set(rh,'visible','off');
      plot(wc(1:j),dc(1:j),'k-','linewidth',2);
      rh = plot(wc(j),dc(j),'r.','markersize',20);    
      FAIRpause(1/100)
    end;
    fprintf('.'); if ~rem(j,50) || j==length(wc), fprintf('\n'); end;
  end;
  fprintf('============== %s done ===================\n',DM{k})
end;
% ==============================================================================