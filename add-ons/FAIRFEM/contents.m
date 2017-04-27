%==============================================================================
% This code is part of the Finite Element Method app for the Matlab-based toolbox
%  FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR/FAIRFEM 
%==============================================================================
%
% Image Registration using Finite Element Method discretizations
%
% This toolbox requires FAIR (obtained from https://github.com/C4IR/FAIRFEM)
%
% Contents of FAIRFEM toolbox
%
% ArmijoDiffeomorghicFEM.m  - backtracked line search, see also ArmijoDiffeomorphic.m
% FAIRplotsFEM.m            - adjusted plot tool
% FEMobjFctn.m              - objective function
% FEMPIRobjFctn.m           - objective function for mass-preserving registration
% MLIRFEM.m                 - multi-level strategy (uniform refinement)
% elasticFEM.m              - elastic regularizer
% hyperElasticFEM.m         - hyperelastic regularizer
% getGradientMatrixFEM.m    - computes gradient operator 
% volTetraGrid.m            - computes volume of elements after deformation
% getDeterminant.m          - computes Jacobian determinant
%
% Meshes: A selection of triangular and tetrahedral meshes can be found in
% the folder meshes
%
% Examples:
%
% EFEM_CompareTriangulations2D.m  - compares different triangulations in 2D
% EFEM_Hands2DMLIRFEM.m           - multi-level registration using FEM
% EFEM_Nodal_vs_FEM2D.m           - compares hyperelastic regularizer
% EFEM_FEMPIRE_Gaussian2D.m       - mass-preserving registration of Gaussians
% EFEM_SSDvsMP.m                  - compares mass-preserving and SSD registration
%
% ==================================================================================

function debit = contents
if nargout == 0, help(mfilename); return; end;

debit = {
          'ArmijoDiffeomorphicFEM.m'
          'FAIRplotsFEM.m'
          'FEMPIREobjFctn.m'
          'FEMobjFctn.m'
          'MLIRFEM.m'
          'contents.m'
          'elasticFEM.m'
          'getDeterminant.m'
          'getGradientMatrixFEM.m'
          'hyperElasticFEM.m'
          'plotTriMesh.m'
          'volTetraGrid.m'
        };
