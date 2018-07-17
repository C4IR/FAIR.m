%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% This files simplifies the generation of FAIR data
% Input
%    example        the name of the game (for figures etc)
%    fileT        name of a image file for the template  image
%    fileR        name of a image file for the reference image
%    varargin    list of additional variables such as omega, m, LM
%==============================================================================

function expfile = jpgs2data(example,fileT,fileR,varargin)

if nargin == 0,
  help(mfilename)
  runMinimalExample
  return;
end

omega     = [];
m         = [];
LM        = [];
viewPara  = {'viewImage','viewImage2D','colormap','gray(256)'};
imgPara   = {'imgModel','linearInter'};
overwrite = 0;

for k=1:2:length(varargin), % overwrite defaults
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

if isempty(example),
  caller  = dbstack;
  expfile = caller(min(2,length(caller))).name;
  example = expfile(6:end-5);
else
  expfile = sprintf('setup%sData',example);
end
expfile = fullfile(FAIRpath,'temp',expfile);

FAIRmessage(expfile)
if exist([expfile,'.mat'],'file') & ~overwrite,
  fprintf('[%s] %s exists\n',mfilename,expfile);
  return
end;

% do whatever needed to be done to get your data here
image = @(str) double(flipud(imread(str))'); % reads and converts

% load the original data, set domain, initial discretization, and grid
dataT = image(fileT);
dataR = image(fileR);

if isempty(omega), omega = [0,size(dataT,1),0,size(dataT,2)];  end;
if isempty(m),     m     = size(dataR);                        end;

viewImage('reset',viewPara{:});
imgModel('reset',imgPara{:});

omegaI = @(i) omega(min(i,size(omega,1)),:);
xc     = @(k) getCellCenteredGrid(omegaI(k),m);

viewData  = @(I,k) viewImage(imgModel(I,omegaI(k),xc(k)),omegaI(k),m);

FAIRfigure(1,'figname',expfile); clf;
subplot(1,2,1); viewData(dataT,1); title('template');
subplot(1,2,2); viewData(dataR,2); title('reference');

ML = getMultilevel({dataT,dataR},omega,m,'fig',2);

% save to outfile
save(expfile,'example','dataT','dataR','omega','m','ML','viewPara','imgPara');
if ~isempty(LM),
  save(expfile,'-append','LM','example');
end;

%------------------------------------------------------------------------------

function runMinimalExample
setup2DhandData;

%==============================================================================

