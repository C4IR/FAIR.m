%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% KERNEL/LANDMARKS
%
%  Contents of FAIR's landmark registration
%
% For an extended documentation, see:
% Jan Modersitzki. FAIR: Flexible Algorithms for Image Registration, SIAM, 2009.
% http://www.siam.org/books/fa06/
% 
% Contents of FAIR landmark based registration Toolbox
%
% contents           - this file
% LMreg                   - the landmark registration toolit, see E5_Hands_TPS for an example
% YTPS                 - computes the landmark based transformation Y
% evalTPS                   - evaluates the Thin-Plate-Spline solution 
%                      based on coefficients
% getLandmarks       - enables the determination of landmarks in 2D
% getTPScoefficients - computes the coefficients for the TPS solution
% plotLM             - convenient way for plotting landmarks
% TPSmex.m           - wrapper for cpp version pf YTPS
%
%  see also E5_Hands_TPS
%
%==============================================================================

function debit = contents
if nargout == 0, help(mfilename); return; end;

debit = {
  'contents.m'
  'LMreg.m'
  'YTPS.m'
  'evalTPS.m'
  'getLandmarks.m'
  'getTPScoefficients.m'
  'plotLM.m'
  'TPSmex.m'
  'TPSmexC.cpp'
  'testLandmarks.m'
  };
%==============================================================================
