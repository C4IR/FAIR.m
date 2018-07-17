%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: Used to create Fig. 3.13 p.40
%
%==============================================================================

setup2DUSData;
omega = ML{end}.omega; m = ML{end}.m;
distance('reset','distance','SSD');
imgModel('reset','imgModel','splineInter');
lMax = length(ML);
m  = @(l) ML{l}.m;
xc = @(l) getCellCenteredGrid(omega,m(l));

Name  = @(c,l) sprintf('ML%c-%d',c,l);
Write = @(c,l,R,m) imwrite(uint8(flipud(reshape(R,m)')),...
  fullfile(FAIRpath,'temp',[Name(c,l),'.jpg']));

for level=3:length(ML);
  R  = imgModel('coefficients',ML{level}.T,[],omega);
  Rc = imgModel(R,omega,xc(level));
  Rf = imgModel(R,omega,xc(lMax));
  
  FAIRfigure(level); clf;
  subplot(1,2,1); viewImage(Rf,omega,m(lMax));
  subplot(1,2,2); viewImage(Rc,omega,m(level)); FAIRpause(1/6);
  Write('f',level,Rf,m(lMax));
  Write('c',level,Rc,m(level));
end;
%==============================================================================
