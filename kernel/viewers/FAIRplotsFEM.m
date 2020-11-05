%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% 
% function FAIRplots(task,varargin);
% 
% FAIRplots: generates FAIR plots, see E6_FAIRplots for an example
%
% Create the following 2-by-3 plot:
% _________________________________________________________________
% | R0                 | T0                 |  Tk                 |
% | titleR0            | titleT0            |  title(Tk)          |
% | showImage(R0)      | showImage(T0)      |  showImage(Tk)      |
% |                    |                    |                     |
% |_______________________________________________________________|
% |                    | D0                 |  Dk                 |
% | titleGrid          | titleT0            |  title(Tk)          |
% | showImage(T0) G0   | showDifference(D0) |  showDifference(Tk) |
% | and showGrid  Gk   |                    |                     |
% |_______________________________________________________________|
%
% The function uses a persitent struct plotOptn controling its behaviour
% and perform the following tasks:
%
%------------------------------------------------------------------------------------
% - FAIRplots('clear')
%   clears all persistent variables and returns
%------------------------------------------------------------------------------------
% - FAIRplots('reset',varargin), FAIRplots('set',varargin)
%   initializes various variables and plot-handles (see details below)
%   overwrites defaults by the varargin list
%   note: no plots are shown using set
%   Example: FAIRplots('reset','mode','NPIR','fig',level);
%
%------------------------------------------------------------------------------------
% - FAIRplots('init',struct('Tc',dataT,'Rc',dataR,'omega',omega,'m',m));
%   initializes the figure 'fig' at position 'position, 
%     R(xc)            -> (2,3,1) title: Rname
%     T(xc)            -> (2,3,2) title: 'T(xc)'
%     T(xc)            -> (2,3,4) title: 'T(xc)'
%     T(xc)-R(xc)      -> (2,3,5) title: '|T(xc)-R|'
%------------------------------------------------------------------------------------
%   para  = struct('Tc',Tstop,'Rc',Rc,'omega',omega,'m',m,'yc',yStop,'Jc',Dc(Tstop));
% - FAIRplots('stop',para); T = Tstop 
%     T(yStop)         -> (2,3,3) title: 'Tstop'
%     T(yStop)-R(xc)   -> (2,3,6) title: '|Tstop-R|=100%'
%     yStop            -> (2,3,4) (update grid)
%------------------------------------------------------------------------------------
%  para.Tc = T0; para.yc = Y0; para.Jstop = para.Jc; para.Jc = Dc(T0);
% - FAIRplots('start',para); 
%     T(y0)            -> (2,3,2) title: Tname
%     T(y0)-R(xc)      -> (2,3,5) title: Dname
%     y0               -> (2,3,4) (update grid)
%------------------------------------------------------------------------------------
%   para.Tc = Tk; para.yc = yk; para.Jc = Dc(Tk); 
%   para.normdY = norm(yk-Y0)/norm(Ystop);
% - FAIRplots(k,para); (note: k is an integer)
%   T(yc)            -> (2,3,3) title: Tanme
%   T(yc)-R(xc)      -> (2,3,6) title: Dname
%   yc               -> (2,3,4) (replace)
%
%------------------------------------------------------------------------------------
%
%  Details on initialization (default values in []):
%
%   the figure is controlled by FIG[=[]] (handle), position[=screen dependent] (positioning),
%   subplots are controllable via handles [TRG]?HANDLE
%   PLOTS[=1] is used to disable the functionality (set plots=0)
%   MODE[=name of calling function] is used to control the default figurename and reference title:
%     mode==PIR*:   
%       figname: PIR*: inter/distance/trafo, dimension/size
%       Rname:   R, m=value of m, length(w)=number of parameters
%     mode==NPIR*:   
%       figname: NPIR*: inter/distance/regularizer, dimension/size
%       Rname:   R, m=value of m, alpha==value of alpha
%     otherwise
%       figname: mode, dimension/size
%       Rname:   R, R, m=value of m
%
%    omega and m are used at various places (computation and output)
%
% The final plots and their entiteling  are performed using reconfigurable function
%  
%    Function(arguments)  default
%    Tshow(T,omega,m)     @viewImage
%    Tname(k)             @(iter) sprintf('T(%d)',iter)
%                         
%    Rshow(R,omega,m)     @viewImage
%    Rname(m)             context dependent, 
%                         e.g. @(m)sprintf('R, %s, length(w)=%d',dimstr(m),wLen);
%                         
%    Dshow(T,R,omega,m)   @(T,R,omega,m) viewImage(255-abs(T-R),omega,m)
%    Dname(k,Jc,Jstop)    @(k,Jc,Jstop) ...
%                         sprintf('|J(%d)/Jstop|=%s%%',k,num2str(abs(100*Jc/Jstop)))
%    
%    Gshow(yc,omega,m)    context dependent, e.g. 
%                         @(yc,omega,m)  plotGrid(yc,omega,m,'spacing',ceil(m/32))
%    Gname(normdY)        @(normdY) sprintf('T(xc), |dY|= %s',num2str(normdY)) 
%    
% run('FAIRplots') or  - see also E6_FAIRplots for an example.
%==============================================================================

%
% FEM Version
%
%
function FAIRplots(task,varargin)

if nargin==0
    help(mfilename)
    runMinimalExample; 
    return;
end

persistent plotOptn
task  = lower(task);

Error  = @(str) error(sprintf('[%s-task=%s] %s\n',mfilename,task,str));
dimstr = @(m) sprintf('[%s]',sprintf(' %d',m));

% ----- CLEAR plotOptn ----------------------------------------------
if strcmp(task,'clear') 
  if nargin == 1,
    clear plotOptions;
    disp('cleared plotOptn')
  else
    Error('nargin>0')
  end;
  return;
end;

% ----- INITIALIZE --------------------------------------------------
if strcmp(task,'reset') | strcmp(task,'set')
  % RESET plotOptn
  if strcmp(task,'reset'),
    plotOptn = [];
  end;

  % INITIALIZE plotOptn
  fields = {
    'fig'
    'plots'
    'figname'
    'position'
    'mode'
    'omega'
    'm'
    'T0handle'
    'Tkhandle'
    'Tshow'
    'Tname'
    'R0handle'
    'Rshow'
    'Rname'
    'D0handle'
    'Dkhandle'
    'Dshow'
    'Dname'
    'G0handle'
    'Gkhandle'
    'Gshow'
    'Gname' 
    };

  for j=1:length(fields)
    if ~isfield(plotOptn,fields{j}),
      plotOptn = setfield(plotOptn,fields{j},[]);
    end;
  end;

  % prepare the figname
  J = find(strcmp(varargin,'mode'));
  if isempty(J),
    mode = dbstack;               
    mode = mode(min(length(mode),2)).name;
  else
    mode = varargin{min(J)+1};
  end;

  % set defaults
  defaults = {
    'mode',     mode,...
    'plots',    1,...
    'position', [],...
    'Rshow',    @viewImage,...
    'Tshow',    @viewImage,...
    'Tname',    @(iter) sprintf('T(%d)',iter),...
    'Dshow',    @(T,R,omega,m) viewImage(255-abs(T-R),omega,m),...
    'Dname',    @(j,Jc,Jstop) ...   
        sprintf('|J(%d)/Jstop|=%s%%',j,num2str(abs(100*Jc/Jstop))),...
    'Gshow',    @(yc,omega,m)  plotGrid(yc,omega,m,'spacing',ceil(m/32)),...
    'Gname',    @(normdY) sprintf('T(xc), |dY|= %s',num2str(normdY)) 
    };
 
% percentage of ratio Jc to Jstop
%     'Dname',    @(j,Jc,Jstop) ... % percentage of ratio difference Jstop and Jc to Jstop
%         sprintf('|J(%d)/Jstop|=%s%%',j,num2str(abs(100*(Jstop-Jc)/abs(Jstop)))),...

  for j=1:2:length(defaults),
    if ~isfield(plotOptn,defaults{j}),
      plotOptn.(defaults{j}) = defaults{j+1};
    elseif isempty(plotOptn.(defaults{j})),
      plotOptn.(defaults{j}) = defaults{j+1};
    end;
  end;

  % collect information about current status of the toolbox
  try, wLen = length(trafo('w0')); catch, wLen        = nan;                end;
  try, imgModelStr = imgModel;     catch, imgModelStr = 'imgModel:void';    end;
  try, trafoStr = trafo;           catch, trafoStr    = 'trafo:void';       end;
  try, distStr  = distance;        catch, distStr     = 'distance:void';    end;
  try, regStr   = regularizer;     catch, regStr      = 'regularizer:void'; end;

  % if necessary TUNE figname
  if ~any(strcmp(varargin,'figname')),
  
    if length(plotOptn.mode) > 2 & strcmp(lower(plotOptn.mode(1:3)),'pir'),
      str = sprintf('%s: %s/%s/%s, ',plotOptn.mode,imgModelStr,distStr,trafoStr);
    elseif length(plotOptn.mode) > 3 & strcmp(lower(plotOptn.mode(1:4)),'npir'),
      str = sprintf('%s: %s/%s/%s, ',plotOptn.mode,imgModelStr,distStr,regStr);
    else
      str = sprintf('%s: ',mode);
    end;
    plotOptn.figname = ...
      @(omega,m) sprintf('%s%dD, m=%s',str,length(omega)/2,dimstr(m));;
  end;

  % if necessary TUNE Rname
  if ~any(strcmp(varargin,'Rname')),    
    if length(plotOptn.mode) > 2 & strcmp(lower(plotOptn.mode(1:3)),'pir'),
      plotOptn.Rname = @(m)sprintf('R, %s, length(w)=%d',dimstr(m),wLen);
    elseif length(plotOptn.mode) > 3 & strcmp(lower(plotOptn.mode(1:4)),'npir'),
      plotOptn.Rname = @(m) sprintf('R, %s, \\alpha=%s',...
        dimstr(m),num2str(regularizer('get','alpha')));
    else
      plotOptn.Rname = @(m) sprintf('R, %s',dimstr(m));
    end;
  end;
  
  % ----- OVERWRITE DEFAULTS ------------------------------------------
  for k=1:2:length(varargin), % overwrites defaults
    if ~isfield(plotOptn,varargin{k}),
      Error(sprintf('field %s not in use',varargin{k}));
    end;
    plotOptn.(varargin{k}) = varargin{k+1};
  end;
  % -------------------------------------------------------------------
  return;
end;
% -------------------------------------------------------------------


if ~plotOptn.plots, return; end;
pause = @(l) builtin('pause',l);
fig = plotOptn.fig;

% extract the parameters for plots
if nargin>1, 
  T      = getField(varargin{1},'Tc'); 
  R      = getField(varargin{1},'Rc'); 
  omega  = getField(varargin{1},'omega'); 
  m      = getField(varargin{1},'m'); 
  yc     = getField(varargin{1},'yc');
  normdY = getField(varargin{1},'normdY');
  Jc     = getField(varargin{1},'Jc');
end;


switch task,

  case 'clear', return; 
  case 'reset', return; 
  case 'set',   return; 
    
  case 'init',
    % activate figure, set colordef and figname
    fig = FAIRfigure(fig);
%     if isempty(fig), fig = figureh; else fig=figureh(fig); end;
    if ~isnumeric(fig),
        fig = fig.Number;
    end;

    plotOptn.fig = fig;
    % extract variables
    [xc,elmat] = getFEMtriangleMesh(omega,m);
    if any(size(T)==1),
      T = reshape(T,m);
      R = reshape(R,m);
    end;
    T  = imgModel(T,omega,xc);
    R  = imgModel(R,omega,xc);

    if ~isempty(plotOptn.position),
      pos = plotOptn.position;
    else
      pos = FAIRposition('fig',fig);
    end;
    figure(fig)
    FAIRfigure(fig,...
      'figname',plotOptn.figname(omega,m),...
      'position',pos);

    if size(omega,2) == 6, % disable Gshow
      plotOptn.Gshow = @(yc,omega,m) [];
    end;

    % plot
    % _____________________________
    % | R       | T     |         |
    % |_________|_______|_________|
    % | T +grid | T -R  |         |
    % |_________|_______|_________|
    figureh(fig); clf
      subplot(2,3,1); 
        plotOptn.R0handle = plotOptn.Rshow(R,xc,elmat);      
        title(plotOptn.Rname(m));
      subplot(2,3,2); 
        plotOptn.T0handle = plotOptn.Tshow(T,xc,elmat);      
        title('T(xc)');
      subplot(2,3,4); 
        plotOptn.G0handle = plotOptn.Tshow(T,xc,elmat);      
        plotOptn.Gkhandle = [];      
        title('T(xc)&grid'); hold on;
        axis manual
      subplot(2,3,5); 
        plotOptn.D0handle = plotOptn.Dshow(T,R,xc,elmat);    
        plotOptn.Dkhandle = [];      
        title('|T(xc)-R|');
        FAIRpause(1/100);
        drawnow;

  case 'stop',
    % plot
    % _____________________________
    % |         |       | Tstop   |
    % |_________|_______|_________|
    % |   +grid |       | Tstop-R |
    % |_________|_______|_________|
    figureh(fig);
      subplot(2,3,3); 
        plotOptn.Tkhandle = plotOptn.Tshow(T,xc,elmat);     
        title('T^{stop}');
      subplot(2,3,6); 
        plotOptn.Dkhandle = plotOptn.Dshow(T,R,xc,elmat);
        title('|T^{stop}-R|, J^{stop}=100%');
      subplot(2,3,4); 
        set(plotOptn.Gkhandle,'visible','off');
        plotOptn.Gkhandle = plotOptn.Gshow(yc,xc,elmat);
    pause(1/100); drawnow;
    plotOptn.Jstop = Jc;
    
  case 'start',
    % plot
    % _____________________________
    % |         | T0    |         |
    % |_________|_______|_________|
    % |   +grid | T0-R  |         |
    % |_________|_______|_________|
    figureh(fig);
      subplot(2,3,2); 
        plotOptn.T0handle = plotOptn.Tshow(T,xc,elmat);    
        title(plotOptn.Tname(0));
      subplot(2,3,5); 
        plotOptn.D0handle = plotOptn.Dshow(T,R,xc,elmat);   
        title(plotOptn.Dname(0,Jc,plotOptn.Jstop));
     subplot(2,3,4); 
       set(plotOptn.Gkhandle,'visible','off');
       plotOptn.Gkhandle = plotOptn.Gshow(yc,xc,elmat);
    pause(1/100); drawnow;
    plotOptn.J0 = Jc;
  otherwise, 
    if ~isnumeric(task),
      warning(['don''t no how to deal task <',task,'>!']);
      return;
    end;
    
    % plot
    % _____________________________
    % |         |       | Tc      |
    % |_________|_______|_________|
    % |   +grid |       | Tc-R    |
    % |_________|_______|_________|    
    figureh(fig)
      subplot(2,3,4); 
        set(plotOptn.Gkhandle,'visible','off');
        plotOptn.Gkhandle = plotOptn.Gshow(yc,xc,elmat); 
        title(plotOptn.Gname(normdY));
      subplot(2,3,3); 
        plotOptn.Tkhandle = plotOptn.Tshow(T,xc,elmat);    
        title(plotOptn.Tname(task));
      subplot(2,3,6); 
        plotOptn.Dkhandle = plotOptn.Dshow(T,R,xc,elmat);      
        title(plotOptn.Dname(task,Jc,plotOptn.Jstop));
    drawnow;    
end;

drawnow
%------------------------------------------------------------------------------

function v = getField(s,field);
if isfield(s,field), v = s.(field); else v = []; end;

%------------------------------------------------------------------------------

function s = setDefault(s,varargin);
for j=1:2:length(varargin),
  if ~isfield(s,varargin{j}),
    s.(varargin{j}) = varargin{j+1};
  elseif isempty(s.(varargin{j})),
    s.(varargin{j}) = varargin{j+1};
  end;
end;

%------------------------------------------------------------------------------

function runMinimalExample
setup2DhandData
imgModel('reset','imgModel','splineInter','regularizer','moments','theta',1e0);
level = 4; omega = ML{level}.omega; m = ML{level}.m;
[T,R] = imgModel('coefficients',ML{level}.T,ML{level}.R,omega(1,:),'out',0);
x0  = getCellCenteredGrid(omega(1,:),m);
x1  = rotation2D( 25*pi/180,x0,'c',(omega(2:2:end)-omega(1:2:end))'/2);
x2  = rotation2D(-25*pi/180,x0,'c',(omega(2:2:end)-omega(1:2:end))'/2);  
Rc  = imgModel(R,omega(end,:),x0);
T0  = imgModel(T,omega(end,:),x0);
T1  = imgModel(T,omega(end,:),x1);
T2  = imgModel(T,omega(end,:),x2);

FAIRfigure(2);
subplot(1,2,1); viewImage2Dsc(T0,xc,elmat,'title','template','colormap','gray')
subplot(1,2,2); viewImage2Dsc(Rc,xc,elmat,'title','reference','colormap','gray')

FAIRplots('reset','mode','testing-PIR','fig',10,'plots',1);
FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m));

para = struct('Tc',T0,'Rc',Rc,'omega',omega,'m',m,'yc',x0,'Jc',100);
FAIRplots('stop',para); 
pause(3);

para = struct('Tc',T1,'Rc',Rc,'omega',omega,'m',m,'yc',x1,'Jc',14.3);
FAIRplots('start',para); 
pause(3);

para = struct('Tc',T2,'Rc',Rc,'omega',omega,'m',m,'yc',x2,'Jc',0.3);
FAIRplots(14,para);
%==============================================================================