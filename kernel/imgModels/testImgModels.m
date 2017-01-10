%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% This is a testing environment for the folder kernel/imgModels
% 1. Based on data/contents, a list of required files is generated and it 
%    is verified, that all files are present; additional files are listed.
% 2. All c-files are compiled.
% 3. All files are executed.
% 4. Check administrative modul imgModel.m
%==============================================================================

folder = fileparts(which(mfilename));
files  = testStart(folder);

%% check the administrative module

%% Test module imgModel, start with syntax
imgModel('reset','imgModel','linearInterMatlab','regularizer','none','theta',1);
imgModel('disp');
[scheme,parameter] = imgModel;
imgModel('clear');

testEnd;
%==============================================================================
