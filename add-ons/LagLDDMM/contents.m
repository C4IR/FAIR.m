%==============================================================================
% This code is part of the Matlab-based toolbox 
% LagLDDDM - A Lagrangian Gauss--Newton--Krylov Solver for Mass- and 
%                        Intensity-Preserving Diffeomorphic Image Registration
% 
% For details and license info see 
% - https://github.com/C4IR/FAIR.m/tree/master/add-ons/LagLDDMM
%==============================================================================
%
% For a quick example, run: EMPLDDMM_3Dmouse_mfDiffusionCC.m
%
% Contents of LagLDDMM toolbox:
%
% GaussNewtonLDDMM.m                            - tailored GaussNewton.m
% LDDMMobjFctn.m                                - intensity-preserving objective function
% MPLDDMMobjFctn.m                              - mass-preserving objective function
% MLLDDMM.m                                     - tailored multilevel
% dctn.m                                        - N-D dct
% eigLaplacian.m                                - computes eigenvalues of diff ops
% linearInterGrid.m                             - grid interpolation
% getLinearInterGridMatrix.m                    - derivative of grid interpolation
% getLinearInterMatrix.m                        - derivative of interpolation
% getPICMatrixAnalyticIntegral.m                - push-forward matrix
% getTrafoFromInstationaryVelocityRK4.m         - integrator for characteristics
% getTrafoFromVelocityRK4.m                     - integrator for characteristics
% getVelocityStartingGuess.m                    - constructs starting guess of correct size
% idctn.m                                       - N-D idct
% mfCurvatureST.m                               - matrix-free spatio-temporal curvature regularizer
% mfDiffusionCC.m                               - matrix-free spatial diffusion regularizer 
% mfDiffusionST.m                               - matrix-free spatio-temporal diffusion regularizer
% mbCurvatureST.m                               - matrix-based spatio-temporal curvature regularizer
% mbDiffusionCC.m                               - matrix-based spatial diffusion regularizer 
% mbDiffusionST.m                               - matrix-based spatio-temporal diffusion regularizer
% solveSpectral.m                               - spectral solver for diff ops
% spectralPrecondPCG.m                          - preconditioner
% 
%
% Examples:
%
% ELDDMM_2Ddisc2C_mfCurvatureST.m
% ELDDMM_2Ddisc2C_mfDiffusionCC.m
% ELDDMM_2Ddisc2C_mfDiffusionST.m
% ELDDMM_2Ddisc2C_mbCurvatureST.m
% ELDDMM_2Ddisc2C_mbDiffusionCC.m
% ELDDMM_2Ddisc2C_mbDiffusionST.m
% ELDDMM_2Dhands_mfDiffusionCC.m
% ELDDMM_Precond_diffusionCC.m
% ELDDMM_Precond_diffusionST.m
% EMPLDDMM_2DGaussian_mfDiffusionCC.m
% EMPLDDMM_2DGaussian_mfDiffusionST.m
% EMPLDDMM_3Dmouse_mfDiffusionCC.m
% EMPLDDMM_3Dmouse_mfDiffusionST.m
% Ex_HyperElasticNIREP.m
% Ex_LDDMMNIREP.m
%
% ==================================================================================

function debit = contents
if nargout == 0, help(mfilename); return; end;

debit = {
          'GaussNewtonLDDMM.m'
		  'LDDMMobjFctn.m'
		  'MLLDDMM.m'
		  'MPLDDMMobjFctn.m'
		  'README.md'
		  'contents.m'
		  'dctn.m'
		  'eigLaplacian.m'
		  'examples'
		  'getLinearInterGridMatrix.m'
		  'getLinearInterMatrix.m'
		  'getPICMatrixAnalyticIntegral.m'
		  'getTrafoFromInstationaryVelocityRK4.m'
		  'getTrafoFromVelocityRK4.m'
		  'getVelocityStartingGuess.m'
		  'idctn.m'
		  'linearInterGrid.m'
		  'mfCurvatureST.m'
		  'mfDiffusionCC.m'
		  'mfDiffusionST.m'
		  'mbCurvatureST.m'
		  'mbDiffusionCC.m'
		  'mbDiffusionST.m'
		  'solveSpectral.m'
		  'spectralPrecondPCG.m'
          'testLagLDDMM.m'
        };
