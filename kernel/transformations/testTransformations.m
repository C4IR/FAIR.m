%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% This is a testing environment for the files in the folder kernel/data
% 1. Based on data/contents, a list of required files is generated and it 
%    is verified, that all files are present; additional files are listed.
% 2. All c-files are compiled.
% 3. All files are executed.
% 4. Check administrative modul Trafo.m
% 5. Transform certain 2D and 3D grids  
%==============================================================================

FAIRcheckFiles(mfilename);

%% Test module trafo, start with syntax
trafo('reset','trafo','splineTransformation2D',...
  'p',[4,5],'omega',[1,1],'m',[6,40]);
trafo('disp');
[scheme,parameter] = trafo;
trafo('clear');
%% test 2D cases

omega  = [0,1,0,2]; 
m      = [6,7]; 
xc     = getCellCenteredGrid(omega,m);
center = (omega(2:2:end)+omega(1:2:end))'/2;

% initialize various transformations
trafos = {
  'affine2D',{},...
  'affine2Dsparse',{},...
  'rigid2D',{},...
  'rotation2D',{'c',center},...
  'translation2D',{},...
  'splineTransformation2D',{'p',[4,5],'omega',omega,'m',m},...
  'splineTransformation2Dsparse',{'p',[4,5],'omega',omega,'m',m}
  };

for k=2:length(trafos)/2,
  fprintf(2,'\n\n\ntest [ %s ] --- %d of %d\n',...
    trafos{2*k-1},k,length(trafos)/2)
  optn = trafos{2*k};
  trafo('reset','trafo',trafos{2*k-1},'debug','on',optn{:});
  trafo('disp');
  wc = trafo('w0');
  fctn = @(wc) trafo(wc,xc);
  checkDerivative(fctn,wc+randn(size(wc)),'fig',1+mod(k,2));
  title(trafo)
  builtin('pause',2);
end;

%% test 3D cases
omega = [0,1,0,2,0,3]; m = [6,7,8]; xc = getCellCenteredGrid(omega,m);
trafos = {
  'affine3D',{},...
  'affine3Dsparse',{},...
  'rigid3D',{},...
  'splineTransformation3Dsparse',{'p',[4,5,6],'omega',omega,'m',m},...
  'translation3D',{}
  };

for k=1:length(trafos)/2,
  fprintf('\n\n\ntest [ %s ] --- %d of %d\n',...
    trafos{2*k-1},k,length(trafos)/2)
  optn = trafos{2*k};
  trafo('reset','trafo',trafos{2*k-1},'debug','on',optn{:});
  trafo('disp');
  wc = trafo('w0');
  fctn = @(wc) trafo(wc,xc);
  checkDerivative(fctn,wc+randn(size(wc)),'fig',1+mod(k,2));
  title(trafo)
  builtin('pause',2);
end;

testEnd;
%==============================================================================
