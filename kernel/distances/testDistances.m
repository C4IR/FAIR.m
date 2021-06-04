%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% This is a testing environment for the files in the folder
% kernel/distances
% 1. Based on data/contents, a list of required files is generated and it 
%    is verified, that all files are present; additional files are listed.
% 2. All c-files are compiled.
% 3. All files are executed.
% 4. Check administrative modul distance.m
% 5. Test explicitly several specific distances
%==============================================================================

FAIRcheckFiles(mfilename);

%% Test module distance, start with syntax
distance('reset','distance','MIcc',...
  'tol',1e-7,'minT',0,'maxT',256,'nT',60,'minR',0,'maxR',256,'nR',60);
distance('disp');
[scheme,parameter] = distance
distance('clear');
distance('disp');

%% 5. setup test cases and run these 
fprintf(' ... %s\n','generate hand data and setup transformation model:')
setup2DhandData
imgModel('reset','imgModel','splineInter','regularizer','moments','theta',1e0);
level = 4; omega = ML{level}.omega; m = ML{level}.m;
[T,R] = imgModel('coefficients',ML{level}.T,ML{level}.R,omega(1,:),'out',0);
xc  = getCellCenteredGrid(omega(1,:),m);
Rc  = imgModel(R,omega(end,:),xc);

FAIRfigure(2);
subplot(1,2,1); viewImage2Dsc(T,omega,m,'title','template','colormap','gray')
subplot(1,2,2); viewImage2Dsc(R,omega,m,'title','reference','colormap','gray')

%% test several distances
distances = {'SSD','NCC','MIspline','MImex','NGFdot'};
clear fig; l = 0;

for k=1:length(distances);

  distance('reset','distance',distances{k});
  switch distance,
    case {'SSD','NCC'},
    case 'MIcc',
      %setup default parameter
      distance('set','tol',1e-7,...
        'minT',0,'maxT',256,'nT',60,'minR',0,'maxR',256,'nR',40);
      distance('disp');
    case 'NGFdot',
      %setup default parameter
      distance('set','edge',100);
      distance('disp');
  end;

  type = {'function','derivative','Gauss-Newton'};
  for k = 1:length(type),
    fctn = @(x) parseDistance(type{k},T,Rc,omega(end,:),m,x);
    l = l+1;
    fig(l)=checkDerivative(fctn,xc);
    ylabel([distance,'-',type{k}])
    builtin('pause',1);
  end;
end;

testEnd;
%==============================================================================

