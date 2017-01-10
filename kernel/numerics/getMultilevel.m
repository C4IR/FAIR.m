%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% function [ML,minLevel,maxLevel,fig] = getMultilevel(IS,omega,m,varargin)
%
% compute a multi-level representation of an image IS or a list of images 
% IS = {T,R,MT,MR,...}
%
%  Input:
%      IS        an image or a list of images
%      omega     domain specification
%      m         number of discretization points
%      varargin  optional parameter, see below
%
% Output:
%   ML       struct containing multi-level representation
%                   ML{j} =  {T_j,R_j,omega,m_j}
%   minLevel     coarsest level
%   maxLevel     finest   level
%   fig          handle to the graphical output
%
% A continuous representation of an image serves as a starting point. Using a cell 
% centered grid, the data is replace by interpolated values. The representation on 
% a coarser level is obtained by averaging over adjacent cells (see the filters for 
% options). 
%
% ML{level} is a structure containing the image(s), omega, and m,
% where level runs 0:maxLevel, ML{level} = [] for level < minLevel.
% Note that we assume m=2^p, p integer; the original data size is arbitrary.
%==============================================================================

function [ML,minLevel,maxLevel,fig] = getMultilevel(IS,omega,m,varargin)

if nargin == 0,    % show help and provide minimal example
  help(mfilename); 
  runMinimalExample;
  ML = 'endOfMinimalExample';
  return;
end;

if nargin == 1,
  % simply return ML, minLevel, maxLevel
  ML = IS; maxLevel = length(ML); fig = 0;
  for minLevel=maxLevel:-1:1, 
    if minLevel == 1 || isempty(ML{minLevel-1}), break; end; 
  end;
  return;
end;
    
% start the work
fig      = 2;                  % figure number for output
dopause  = 0;                  % make pause for demonstrations
filter   = 'block';            % discrete smoothing kernel
minLevel = 3;                  % minimal level size
names    = {'T','R','Q','Q1','Q2','Q3','Q4','Q5'};
restrictdim = ones(size(m));   % by default: restrict all dimensions
imgModel = @linearInter;       % use linear interpolation as a default
                               % image model 

for k=1:2:length(varargin),    % overwrite defaults  
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

if ~iscell(IS), IS = {IS}; end;% in case of a single image input
lenIS = length(IS);                 % number of images to be handled
dim   = length(omega)/2;            % spacial dimension
maxLevel = ceil(log2(min(m)));      % finest level

% messaging 
msg = sprintf('%s(%d image(s), filter=%s, %dD data,level=%d:%d)',...
  mfilename,lenIS,filter,dim,minLevel,maxLevel);

fprintf('%s, figure=%d\n',msg,fig);
fprintf('note: model from data of any dimension is sampled on a grid of size m\n');
fprintf('      at present it is assume that m = 2.^p, p an interger\n');
fprintf('      ML is an array of cells, where the i-th entry contains the data\n');
fprintf('      of size 2^p, with some modifications is not all dimensions are equal\n');

% some output to the console
fprintf('%s: [',mfilename);

omegak = @(k)omega(min(1+(k>1),size(omega,1)),:);

% start the loop from maxlevel to minlevel
for level = maxLevel:-1:minLevel
  if level == maxLevel,
    % replace data by the interpolant on x
    for k=1:lenIS,
      sizeData = size(IS{k});
      % note: data is either m1-by-....-by-md 
      %                   or m1-by-....-by-md-by-nrVolume
      % but the 1d case is tricky

      if (dim>1) && (length(sizeData) > dim),
        nrVolume = sizeData(end), 
        sizeData(end) = [];
        keyboard
      else
        nrVolume = 1;
      end;
                    
      data  = reshape(IS{k},prod(sizeData),nrVolume);
      block = zeros(prod(m),nrVolume);      
      xc    = getCellCenteredGrid(omegak(k),m);
      % create dsample of size m from image model (interpolation) obtained from data
      for vol=1:size(data,2),
        block(:,vol) = imgModel(reshape(data(:,vol),sizeData),omegak(k),xc);
      end
      
      IS{k} = reshape(block,[m,nrVolume]);
    end;
    fprintf('%d',level);
  else
    % restrict the image to the coarser grid
    L = ML{level+1};
    for k=1:lenIS, % run over all images
      for j=1:dim, % run over all dimensions
        if restrictdim(j),
            IS{k} = restrict(IS{k},dim,j,filter); 
        end
      end; 
    end;
    fprintf(',%d',level);
  end;
  
  % store the current level data in a struct
  L.m = size(IS{1}); L.m = L.m(1:dim); L.omega = omega; 
  for k= 1:lenIS,
    L = setfield(L,names{k},IS{k});
  end;
  ML{level} = L;
  
  if fig, % do some plots
    if level == maxLevel,
       FAIRfigure(fig,'figname',msg,'position','default');
    end;
    
    str = @(k) sprintf('%s(level=%d), %s',...
      names{k},level,sprintf('m=[%s]',sprintf(' %d',L.m)));
    p0 = level-minLevel+1; 
    dp = (maxLevel-minLevel+1);
    figureh(fig);
    for k=1:lenIS,
      xc = getCellCenteredGrid(omegak(k),L.m);
      subplot(lenIS,dp,p0+(k-1)*dp);    
      viewImage(imgModel(IS{k},omegak(k),xc),omegak(k),L.m); 
      title(str(k));
    end;
    if dopause, pause; else drawnow; end;
  end;  
end;
fprintf('] done\n');

%------------------------------------------------------------------------------
function T = restrict(T,dim,j,filter)
J = [j,setdiff(1:dim,j)]; J = [J dim+1]; % bring dimension j to front
T = permute(T,J);
if rem(size(T,1),2), T(end+1,:,:,:) = T(end,:,:,:); end;
switch filter,
  case 'gaussian', 
    Ta = zeros(size(T)+[2,zeros(1,ndims(T)-1)]);
    Ta(2:end-1,:,:,:) = T; Ta([1,end],:,:,:) = Ta([2,end-1],:,:,:); 
    T =   0.125 * Ta(1:end-3,:,:,:) + 0.375 * Ta(2:end-2,:,:,:) ...
        + 0.375 * Ta(3:end-1,:,:,:) + 0.125 * Ta(4:end,:,:,:);
    T = T(1:2:end,:,:,:);
  case 'block'
    T = (T(1:2:end,:,:,:)+T(2:2:end,:,:,:))/2; 
  case 'harmonic'
    T = 1./(T+eps);
    T = (T(1:2:end,:,:,:)+T(2:2:end,:,:,:))/2; 
    T = 1./T;
end;
T = ipermute(T,J); % bring j back home

%------------------------------------------------------------------------------

function runMinimalExample
Tdata = double(flipud(imread('hands-T.jpg'))');
Rdata = double(flipud(imread('hands-R.jpg'))');
omega = [0,20,0,25]; % specify physical domain
m     = size(Tdata);

% set view options and interpolation options and initialize viewer and interpolator
viewPara = {'viewImage','viewImage2D','colormap','gray(256)'};
viewImage('reset',viewPara{:});
% create multilevel representation of the data
ML = getMultilevel({Tdata,Rdata},omega,m,'fig',2);

%==============================================================================
