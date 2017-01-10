%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% KERNEL/REGULARIZERS
%
% Contents of FAIR's Reguarization Toolbox
% 
% general purpose tools:
%   contents                - this file
%   regularizer             - the specific regularizer used in FAIR
%     initialize: regularizer('reset','regularizer','mbElastic',...
%                             'alpha',1e3,'mu',1,'lambda',0);
%     use:        [Sc,dS,d2S] = regularizer(yc-yRef,omega,m);
%
%  The regularizers for NPIR are:
%  ------------------------------
%  curvature                curvature regularizer (cell-centered grid, requires parameter alpha)
%                           see, e.g., E10_2Ddisc2C_curvature
%  curvaturemexC            matrix free C version of curvature regularizer
%  curvatureDiagMex         matrix free access to the diagonal of the operator
%  curvatureHessianMex      matrix free access to the second derivative of the operator
 
%  elastic                  (linear) elastic regularizer (staggered grid)
%                             initialize: 'alpha', 'mu', 'lambda'
%                             see, e.g., E10_2Ddisc2C_elastic
%  hyperelastic             hyperelastic regularizer    (nodal grid)
%                             initialize: 'alphaLength', 'alphaArea', 'alphaVolume'
%                             see, e.g., E10_2Ddisc2C_hyperElastic
%
%  (linear) Differential operators are built in:
%
%  getCurvatureMatrix       generates curvature regularizer matrix (cell-centered grid)
%  getElasticMatrixNodal    generates elastic regularizer matrix (nodal grid)
%  getElasticMatrixStg      generates elastic regularizer matrix (staggeres grid)
%  getGradientNodal         generates gradient operator matrix (nodal grid)
%
%  For (nonlinear) hyperelastic regularization area and volume of tetrahedral partition 
%  are computed in
%  geometry                contains areas and volumes functions
%                          used in hyperElastic
%  geometrymexC.h / .cpp   MEX version of area and volume computation
%
%
% =============================================================================

function debit = contents
if nargout == 0, help(mfilename); return; end;

debit = {
          'contents.m'

          'regularizer.m'
          
          'elastic.m'
          'getElasticMatrixNodal.m'
          'getElasticMatrixStg.m'
          
          'curvature.m'
          'curvatureMexC.cpp'
          'getCurvatureMatrix.m'
          'curvatureDiagMex.cpp'
          'curvatureHessianMex.cpp'       

          'geometry.m'
          'geometryMexC.cpp'
          'geometryMexC.h'
          'geometryC.cpp'
          'geometryC.h'
          'hyperElastic.m'
          
          'getGradientNodal.m'
          'testRegularizers.m'
  };
%==============================================================================
