%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% 
% KERNEL/DISTANCES
%
% Contents of FAIR's Distance Toolbox
%
% This folder provides several distance measures, wrappers and utilities.
%
% general purpose tools:
%   contents        - this file
%   distance        - the specific distance measure used in FAIR
%     initialize:   distance('reset','distance','SSD');
%     usage:        D = distance(Tc,Rc,omega,m);
%   testDistances   - test the files in this folder
%   see also E9_Hands_MLIR_SSD_mbElas
%
% measures and wrappers:
%   MI              - wrapper for MIcc
%   MImex           - Mutual Information based on rhoSplineC (C-implementation)
%   MIspline        - Mutual Information based on rhoSpline  (MATLAB-implementation)
%   rhoSpline       - MATLAB implementation of joint entropy estimator (spline kernel)  
%   rhoSplineC      - C based implementation of joint entropy estimator (spline kernel)  
%   NCC             - Normalized Cross Correlation
%   NCCmex          - Wrapper for NCCmexC
%   NCCmexC         - C version of NCC
%   NGF             - wrapper for NGFdot
%   NGFdot          - Normalized Gradient Fields, dot product based (THE implementation)
%   NGFcross        - Normalized Gradient Fields, cross product based (to be used with care)
%   NGFmex          - Wrapper for NGFdotMexC
%   NGFdotMexC      - C version of NGFdot
%   SSD             - Sum of Squared Differences 
%   SSDmex          - Wrapper for SSDmexC
%   SSDmexC         - C version of SSD
%   parseDistance   - parse distances for derivative checks
%  see also E9_Hands_MLIR_SSD_mbElas
% 
%==============================================================================

function debit = contents
if nargout == 0, 
    help(mfilename); 
    return; 
end;



debit = {
  'contents.m'
  'distance.m'
  'testDistances.m'

  'SSD.m'
  'SSDmex.m'
  'SSDmexC.cpp'
  
  'NCC.m'
  'NCCmex.m'
  'NCCmexC.cpp'
  
  'NGF.m'
  'NGFcross.m'
  'NGFdot.m'
  'NGFmex.m'                                
  'NGFdotMexC.cpp'                           
  
  'rhoSpline.m'
  'rhoSplineC.cpp'
  
  'MI.m'
  'MImex.m'
  'MIspline.m'
  
  'parseDistance.m'
  };
  %==============================================================================
  