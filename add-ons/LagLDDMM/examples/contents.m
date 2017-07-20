%==============================================================================
% This code is part of the Matlab-based toolbox 
% LagLDDDM - A Lagrangian Gauss--Newton--Krylov Solver for Mass- and 
%                        Intensity-Preserving Diffeomorphic Image Registration
% 
% For details and license info see 
% - https://github.com/C4IR/FAIR.m/tree/master/add-ons/LagLDDMM
%==============================================================================
%
% Examples for LagLDDMM toolbox.
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
          'ELDDMM_2Ddisc2C_mfCurvatureST.m'
          'ELDDMM_2Ddisc2C_mfDiffusionCC.m'
          'ELDDMM_2Ddisc2C_mfDiffusionST.m'
		  'ELDDMM_2Ddisc2C_mbCurvatureST.m'
          'ELDDMM_2Ddisc2C_mbDiffusionCC.m'
          'ELDDMM_2Ddisc2C_mbDiffusionST.m'
		  'ELDDMM_2Dhands_mfDiffusionCC.m'
          'ELDDMM_Precond_diffusionCC.m'
          'ELDDMM_Precond_diffusionST.m'
          'EMPLDDMM_2DGaussian_mfDiffusionCC.m'
          'EMPLDDMM_2DGaussian_mfDiffusionST.m'
          'EMPLDDMM_3Dmouse_mfDiffusionCC.m'
          'EMPLDDMM_3Dmouse_mfDiffusionST.m'
          'contents.m'
        };
