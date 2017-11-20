%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% load ONE US image, used as an example for transformations and such
%==============================================================================

checkDataFile
if expfileExists, return; end;

example = '2DUS.jpg';

% do whatever needed to be done to get your data here
image = @(str) double(flipud(imread(str))); % reads and converts

% load the original data, set domain, initial discretization, and grid
dataT = image('US.jpg');
dataR = '';
omega = [0,size(dataT,1),0,size(dataT,2)];
m     = 128*[3,2];
% set view options
viewPara = {'viewImage','viewImage2D','colormap','bone(256)'};
viewImage('reset',viewPara{:});

% set interpolation options
imgPara = {'imgModel','linearInter'};
imgModel('reset',imgPara{:});

xc = getCellCenteredGrid(omega,m);
Tc = imgModel(dataT,omega,xc);
FAIRfigure(1,'figname',mfilename); clf;
viewImage(Tc,omega,m);
title(example,'interpreter','none');

ML = getMultilevel(dataT,omega,m,'fig',2);

% save to outfile
save(expfile,'dataT','dataR','omega','m','ML','viewPara','imgPara');
%==============================================================================
