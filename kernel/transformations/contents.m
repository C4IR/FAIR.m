%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% 
% KERNEL/TRANSFORMATIONS
%
% Contents of FAIR's Transformation Toolbox
%
% contents                     - this file
% trafo                        - the specific transformation in FAIR
%	
%  initialize:   trafo('reset','trafo','rigid2D');
%  usage:        w0      = trafo('w0'); % returns parameterization of identity
%                [yc,dy] = trafo(wc,xc;
%
% testOneTransformation        - tests one transformation model
% testTransformations          - tests the files in this folder
%
% affine2D                     - affine linear transformation for 2D
% affine2Dsparse               - affine linear transformation for 2D (memory efficient)
% affine3D                     - affine linear transformation for 3D
% affine3Dsparse               - affine linear transformation for 3D (memory efficient)
% matVecQw                     - matrix-vector multiplication for efficient versions
% rigid2D                      - rigid  transformation for 2D
% rigid2Dsparse                - rigid  transformation for 2D (memory efficient)
% rigid3D                      - rigid  transformation for 3D
% rigid3Dsparse                - rigid  transformation for 3D (memory efficient)
% rotation2D                   - rotation for 2D
% splineTransformation2D       - spline transformation for 2D
% splineTransformation2Dsparse - spline transformation for 2D (memory efficient)
% splineTransformation3Dsparse - spline transformation for 3D (memory efficient)
% translation2D                - translation for 2D
% translation3D                - translation for 3D
% tensorProdC.c                - C implementation of kron(Q3,Q2,Q1) * w
%                                see splineTransformation3Dsparse.m
%  
%  see also E6_HNSP_PIR_SSD_rotation2D_level4 
%==============================================================================

function debit = contents
if nargout == 0, help(mfilename); return; end;

debit = {
  'contents.m'
  'trafo.m'
  'testOneTransformation.m'
  'testTransformations.m'
  
  'affine2D.m'
  'affine2Dsparse.m'
  'affine3D.m'
  'affine3Dsparse.m'
  'matVecQw.m'
  
  'rigid2D.m'
  'rigid2Dsparse.m'
  'rigid3D.m'
  'rigid3Dsparse.m'
  
  'rotation2D.m'
  
  'splineTransformation2D.m'
  'splineTransformation2Dsparse.m'
  'splineTransformation3Dsparse.m'
  
  'translation2D.m'
  'translation3D.m'
  'tensorProdC.c'
  
  };
%==============================================================================
