%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% 
% 
% KERNEL/TOOLS
%
% Contents of FAIR TOOLS
%
%   contents.m          this file
%   dealOptions.m       deals persistent options for use in main modules
%   reportStatus.m      reports the current setting of the toolbox
%
%   testFAIR.m          tests the complete FAIR toolbox
%   testTools.m         test the kernel/tools folder
%               
%   FAIRtestPara.m      handles the variables used for testings
%   testStart.m         initializes the testing environment
%   testEnd.m           finishes    the testing environment
%               
%   FAIRcheckFolder.m   generates list of files to be tested
%   FAIRcheckFiles.m    checks all files from a list
%   FAIRmake.m          convenient way of calling the matlab mex compiler
%   FAIRcell2struct.m   converts cell to struct and vice versa
%   FAIRmessage.m       generates nice message
%   FAIRinput.m         input, disabled in the testing framework
%   FAIRpause.m         pause, disabled in the testing framework
%   FAIRpath.m          shortcut to the root of the toolbox
%==============================================================================

function debit = contents
if nargout == 0, help(mfilename); return; end;

debit = {
  'contents.m'
  'dealOptions.m'
  'reportStatus.m'

  'testFAIR.m'
  'testTools.m'
  
  'FAIRtestPara.m'
  'testStart.m'
  'testEnd.m'
  
  'FAIRcheckFolder.m'
  'FAIRcheckFiles.m'
  'FAIRmake.m'

  'FAIRcell2struct.m'
  'FAIRmessage.m'
  'FAIRinput.m'
  'FAIRpause.m'
  'FAIRpath.m'
  
};
%==============================================================================

