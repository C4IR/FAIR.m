%==============================================================================
% (c) Jan Modersitzki 2011/07/20, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% For an extended documentation, see:
% Jan Modersitzki. FAIR: Flexible Algorithms for Image Registration, SIAM, 2009.
% http://www.siam.org/books/fa06/
% 
% Contents of FAIR INTERPOLATION TOOLBOX
%
% general purpose tools:
%   contents              - this file
%   inter                 - the specific interpolation scheme used in FAIR
%            initialize: inter('reset','inter','linearIOnter1D');
%            use:        [Tc,dT] = inter(T,omega,Y);
%   getSplineCoefficients - computes spline coefficients
%   motherSpline:         - 1D normalized basis spline and its derivative
%
%   schemes are available for dimension dim, where dim=1,2,3,
% 
%   nnInter               - next neighbor interpolation schemes (not recommended for optimization)
%   linearInter           - linear interpolation schemes
%   linearInterMatlab     - wrapper for MATLAB's interp(dim)
%   linearInterSmooth     - linear approximation schemes
%   splineInter           - spline interpolation schemes
%   cubicInter            - fast local cubic spline interpolation schemes
%
%   mex versions of selected schemes:
%   nnInterMex            - wrapper for C-code nnInterMexC (not recommended for optimization)
%   nnInterMexC           - CPP version of next neighbor interpolation (not recommended for optimization)
%   linearInterMex        - wrapper for C-code linearInterMexC
%   linearInterMexC       - CPP version of linear interpolation
%   linearInterMexNeumannBC  - with Neumann boundary conditions
%   linearInterMexNeumannBCC - with Neumann boundary conditions
%   linearInterSmoothMex  - wrapper for C-code linearInterSmoothMexC
%   linearInterSmoothMexC - CPP version of smooth linear approximation
%   splineInterMex        - wrapper for C-code splineInterMexC.c
%   splineInterMexC       - CPP version of spline interpolation
%   cubicInterMex         - wrapper for C-code cubicInterMexC.c
%   cubivInterMexC        - CPP version of local cubic spline interpolation
%
%  see also E9_Hands_MLIR_SSD_mbElas and BigTutorialInter
% version 2015/05/20
%==============================================================================

function debit = contents
if nargout == 0, help(mfilename); return; end;

debit = {
  'contents.m'
  'checkInterpolation.m'
  'inter.m'
  'nnInter.m'
  'nnInterMex.m'
  'nnInterMexC.cpp'
  'linearInter.m'
  'linearInterMex.m'
  'linearInterMexC.cpp'
  'linearInterMexNeumannBC.m'
  'linearInterMexNeumannBCC.cpp'
  'linearInterSmooth.m'
  'linearInterSmoothMex.m'
  'linearInterSmoothMexC.cpp'
  'linearInterMatlab.m'
  'getSplineCoefficients.m'
  'motherSpline.m'
  'splineInter.m'
  'splineInterMex.m'
  'splineInterMexC.cpp'
  'cubicInter.m'                            
  'cubicInterMex.m'                         
  'cubicInterMexC.cpp'                      

  };
