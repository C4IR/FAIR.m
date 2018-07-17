%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% 
% KERNEL/IMGMODEL
%
% Contents of FAIR's IMAGE MODELS 
%
% general purpose tools:
%   contents              - this file
%   imgModel              - the specific image model used in FAIR
%            initialize: imgModel('reset','imgModel','linearInter1D');
%            use:        [Tc,dT] = imgModel(T,omega,Y);
%   getSplineCoefficients - computes spline coefficients
%   motherSpline:         - 1D normalized basis spline and its derivative
%
%   schemes are available for dimension dim, where dim=1,2,3,
% 
%   nnInter               - next neighbor interpolation schemes 
%                           (not recommended for optimization)
%   linearInter           - linear interpolation schemes
%   linearInterMatlab     - wrapper for MATLAB's interp(dim)
%   linearInterSmooth     - linear approximation schemes
%   splineInter           - spline interpolation/approximation schemes
%   cubicInter            - fast local cubic spline interpolation schemes
%
%   mex versions of selected schemes:
%   nnInterMex            - wrapper for C-code nnInterMexC 
%                           (not recommended for optimization)
%   nnInterMexC           - CPP version of next neighbor interpolation 
%                           (not recommended for optimization)
%   linearInterMex        - wrapper for C-code linearInterMexC
%   linearInterMexC       - CPP version of linear interpolation
%   linearInterMexNeumannBC  - with Neumann boundary conditions
%   linearInterMexNeumannBCC - with Neumann boundary conditions
%   linearInterSmoothMex  - wrapper for C-code linearInterSmoothMexC
%   linearInterSmoothMexC - CPP version of smooth linear approximation
%   splineInterMex        - wrapper for C-code splineInterMexC.c
%   splineInterMexC       - CPP version of spline interpolation
%   cubicInterMex         - wrapper for C-code cubicInterMexC.c
%   cubicInterMexC        - CPP version of local cubic spline interpolation
%
%  see also E9_Hands_MLIR_SSD_mbElas
%==============================================================================

function debit = contents
if nargout == 0, help(mfilename); return; end;

debit = {
  'contents.m'
  
  'imgModel.m'
  'nnInter.m'
  'nnInterMex.m'
  'nnInterMexC.cpp'
  
  'linearInter.m'
  'linearInterMex.m'
  'linearInterMexC.cpp'
  
  'linearInterSmooth.m'
  'linearInterSmoothMex.m'
  'linearInterSmoothMexC.cpp'
  
  'linearInterMatlab.m'
  
  'getSplineCoefficients.m'
  'motherSpline.m'
  'splineInter.m'
  'splineInterMex.m'
  'splineInterMexC.cpp'
  
  'testOneImgModel.m'
  'testImgModels.m'
  };
%==============================================================================

