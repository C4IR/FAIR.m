%==============================================================================
% This code is part of the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% KERNEL/DATA
%
% Contents of FAIR's Numeric Toolbox
% 
% General Purpuse Tools
%   contents            - this file
%   checkDerivative     - tests the implementation of a derivative
%   testNumerics        - test this part of the toolbox
%
% Grid operations (typically for a domain specified by omega of m points)
%   getCellCenteredGrid - generates cell-centered grid 
%   getNodalGrid        - generates nodal grid 
%   getStaggeredGrid    - generates staggered grid 
% 
%   center              - converts any grid to cell-center
%   grid2grid           - converts one grid to another
%   stg2center          - operator for grid mapping
%   nodal2center        - operator for grid mapping
%   nodal2centerC.c     - c-implementation for faster evaluation
%   
% % cell-centered          nodal                  stagged-X1               staggered-X2
% +-----+-----+-----+    n-----n-----n-----n    +-----+-----+-----+    +--x--+--x--+--x--+
% |     |     |     |    |     |     |     |    |     |     |     |    |     |     |     |
% |  o  |  o  |  o  |    |     |     |     |    x     x     x     x    |     |     |     |
% |     |     |     |    |     |     |     |    |     |     |     |    |     |     |     |
% +-----+-----+-----+    n-----n-----n-----n    +-----+-----+-----+    +--x--+--x--+--x--+
% |     |     |     |    |     |     |     |    |     |     |     |    |     |     |     |
% |  o  |  o  |  o  |    |     |     |     |    x     x     x     x    |     |     |     |
% |     |     |     |    |     |     |     |    |     |     |     |    |     |     |     |
% +-----+-----+-----+    n-----n-----n-----n    +-----+-----+-----+    +--x--+--x--+--x--+
%  
% Objective functions
%   PIRBFGSobjFctn      - standard Parameteric Image Registration, BFGS style
%   PIRobjFctn          - standard Parameteric Image Registration, Gauss-Newton style
%   NPIRobjFctn-        - standard Non-Parameteric Image Registration, BFGS style
%   NPIRBFGSobjFctn     - standard Non-Parameteric Image Registration, Gauss-Newton style
%
% Optimization Schemes and Tools
%   optPara             - sets default parameters for various optimization schemes such as
%   SteepestDescent     - steppest descent scheme
%   Nesterov            - Nesterov's scheme
%   GaussNewton         - Gaus-Newton scheme
%   TrustRegion         - Trust region scheme
%   lBFGS               - limited memory BFGS scheme
%   Armijo              - line search scheme
%   ArmijoDiffeomorphic - line search scheme with diffeomorphic constraint
%   solveLinearSystem   - simplifies the solution of the linear systems for Gauss-Newton
%   
% Multilevel Tools
%   getMultilevel       - generates a multi-lrevel representation of the data
%   MLIR                - multi-level image registration, the working horse
%   MLPIR               - multi-level parametric image registration
%==============================================================================

function debit = contents
if nargout == 0, help(mfilename); return; end;

debit = {
  'contents.m'
  'checkDerivative.m'
  
  'center.m'
  'getCellCenteredGrid.m'
  'getNodalGrid.m'
  'getStaggeredGrid.m'
  'grid2grid.m'
  'stg2center.m'
  'nodal2center.m'
  'nodal2centerC.cpp'
  
  'PIRBFGSobjFctn.m'
  'PIRobjFctn.m'
  'NPIRobjFctn.m'
  'NPIRBFGSobjFctn.m'
  
  'optPara.m'
  'SteepestDescent.m'
  'Nesterov.m'
  'GaussNewton.m'
  'solveLinearSystem.m'
  'Armijo.m'
  'ArmijoDiffeomorphic.m'
  'TrustRegion.m'
  'lBFGS.m'

  'getMultilevel.m'
  'MLIR.m'
  'MLPIR.m'

  'testNumerics.m'
  };

%=============================================================================='
