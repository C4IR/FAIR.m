%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% tests the image model coded in file
% runs 1D, 2D and 3D example,
% setup image data Tdata, image domain omega, a discretization m
% and some interpolation grid xc
% computes coefficients Tcoef of representation for image model and 
% [Tc,dT] = feval(file,Tcoef,omega,xc);
% visualize results and perform a derivative check
%==============================================================================
function testOneImgModel(file)


if nargin == 0, 
  linearInterMatlab
  return
end;

FAIRmessage(sprintf('%s // %s',mfilename,file),'-');

if ~strcmp(file,'splineInter')
  coeff = @(Tdata,dim) Tdata;
else
  coeff = @(Tdata,dim) getSplineCoefficients(Tdata,'dim',dim);
end

fprintf(2,'%s: test implementation of <%s>, run minimal examples\n',...
  mfilename,file)

fprintf(2,'check <%s> on 1D example\n',file);
Tdata     = [0,1,4,1,0];
omega     = [0,10];
m         = length(Tdata);
Xdata     = getCellCenteredGrid(omega,m);
Tcoef     = coeff(Tdata,1);

xc      = linspace(-1,11,101);
[Tc,dT] = feval(file,Tcoef,omega,xc);

FAIRfigure(1);
subplot(1,2,1);
plot(xc,Tc,'b-',Xdata,Tdata,'ro');
title(sprintf('%s %d-dim',file,1));
subplot(1,2,2);
spy(dT);
title('dT')


fprintf(2,'check <%s> on 2D example\n',mfilename);
omega = [0,10,0,8];
Tdata = [1,2,3,4;1,2,3,4;4,4,4,4]; m = size(Tdata);
Tcoef = coeff(Tdata,2);
Xdata = getCellCenteredGrid(omega,m);
xc    = getCellCenteredGrid(omega+[-1 1 -1 1],5*m);
[Tc,dT] = feval(file,Tcoef,omega,xc);
DD = reshape([Xdata;Tdata(:)],[],3);
Dc = reshape([xc;Tc],[5*m,3]);

FAIRfigure(2); clf;
subplot(1,2,1);  surf(Dc(:,:,1),Dc(:,:,2),Dc(:,:,3));  hold on;
plot3(DD(:,1),DD(:,2),DD(:,3),'r.','markersize',40); hold off;
title(sprintf('%s %d-dim',file,2));
subplot(1,2,2); spy(dT);                     
title('dT')


fprintf(2,'check <%s> on 3D example\n',mfilename);
omega = [0,1,0,2,0,1]; m = [13,16,7];
Xdata    = getCellCenteredGrid(omega,m);
Y     = reshape(Xdata,[m,3]);
Tdata = (Y(:,:,:,1)-0.5).^2 + (Y(:,:,:,2)-0.75).^2 + (Y(:,:,:,3)-0.5).^2 <= 0.15;
Tcoef     = coeff(Tdata,3);
xc    = getCellCenteredGrid(omega,4*m);
[Tc,dT] = feval(file,Tcoef,omega,xc);

FAIRfigure(3); clf;
subplot(1,2,1); imgmontage(Tc,omega,4*m);
title(sprintf('%s %d-dim',file,3));
subplot(1,2,2); spy(dT);                 
title('dT')

fprintf(2,'run derivative test for <%s> (note: needs to fail)\n',mfilename);
fctn = @(xc) feval(file,Tcoef,omega,xc);
xc   = xc + rand(size(xc));
FAIRfigure(4);
checkDerivative(fctn,xc,'fig',4)

FAIRmessage('-')
%==============================================================================

