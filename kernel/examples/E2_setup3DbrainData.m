%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% Template for setup data
%     (here: use brain MRI's; see also setup3DbrainsData.m for references)
%
% This file initializes the following data (which can then be saved): 
%   dataT       template  image, a d-array of size mD, 
%   dataR       reference image, a d-array of size nD, 
%   omega       domain specification
%             omega = (omega(1),omega(2)) x  (omega(3),omega(24))
%   m           initial discretization 
%   ML          multi-level representation of the data
%   LM        landmarks, if available
%
% For representation and visualization
%   viewPara      options for image viewer
%   imgPara     options for image interpolation
% see also setup3DbrainsData and E2_setupHandsData
%==============================================================================

% load 3D data
load brain3D; 
whos
viewPara  = {'viewImage','imgmontage'};  
viewImage('reset',viewPara{:});

imgPara = {'imgModel','linearInter'};   
imgModel('reset',imgPara{:});
FAIRfigure(2); colormap(gray(256));
ML    = getMultilevel({dataT,dataR},omega,m);
save(fullfile(FAIRpath,'temp','brain3DML.mat'),...
	'dataT','dataR','omega','m','ML','viewPara','imgPara');
%==============================================================================

