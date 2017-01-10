%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% setup test scenarios for a transformation and run the model with options
%==============================================================================
function testOneTransformation(model,options)


if nargin == 0,
  help(mfilename)
  testOneTransformation('affine2D',{});
  return;
end;

if findstr(model,'2D')
  omega = [0,1,0,2];
  m = [6,7]
elseif findstr(model,'3D')
  omega = [0,1,0,2,0.1,3];
  m = [6,7,8]
else
  error('unknown dimension')
end;

xc = getCellCenteredGrid(omega,m);
center = (omega(2:2:end) -omega(1:2:end))/2;

trafo('reset','trafo',model,options{:},'center',center);
trafo('disp');
wc = trafo('w0');
fctn = @(wc) trafo(wc,xc)
%==============================================================================
