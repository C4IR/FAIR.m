%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% Artifical data: boxes in 3D
%==============================================================================

checkDataFile
if expfileExists, return; end;

viewPara = {'viewImage','imgmontage','direction','-zyx','colormap','bone(256)'};
[viewer,viewPara] = viewImage('reset',viewPara{:});

omega = [0,1,0,1,0,1];
m     = [64,64,64];
xc    = getCellCenteredGrid(omega,m);
xc    = reshape(xc,[],3);
dataT = ...
  (abs(xc(:,1)-0.5) < 0.25) ...
  & (abs(xc(:,2)-0.5) < 0.20)...
  & (abs(xc(:,3)-0.5) < 0.22);
xc    = xc(:);
dataT = 200*reshape(dataT,m);
wc    = reshape([eye(3),[-0.1;0;0]]',[],1);
yc    = affine3D(wc,xc);
dataR = reshape(linearInter(dataT,omega,yc),m);

FAIRfigure(1); clf;
subplot(1,2,1); viewImage(dataT,omega,m);
subplot(1,2,2); viewImage(dataR,omega,m);

ML = getMultilevel({dataT,dataR},omega,m,'fig',2);

imgPara   = {'imgModel','linearInter'};
traPara   = {'trafo','affine3D'};
disPara = {'distance','SSD'};
regPara = {'regularizer','mfElastic','alpha',500,'mu',1,'lambda',0};

% save to outfile
save(expfile,'dataT','dataR','omega','m',...
  'ML','viewPara','imgPara','traPara','disPara','regPara');
checkDataFile
%==============================================================================
