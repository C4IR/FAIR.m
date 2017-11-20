%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% KERNEL/EXAMPLES
%
% This folder provides tons of examples that come with the toolbox. 
% The numbering is related to chapters of the FAIR book.
%
%------------------------------------------------------------------------------
% GENERAL FILES
%   contents                            this file
%   spline1D                            codes cubic B-spline in 1D, see E3_bspline for an example
%   testExamples                        testing all example files; see testFAIR for more details
%------------------------------------------------------------------------------
% DATO I/O, IMAGE MODELS, INTERPOLATION
%   E2_setup2DhandData                  setup 2D hand data, initializes viewer,
%                                       computes and visualizes multi-level representation
%   E2_setup3DbrainData                 setup 3D brain data, initializes viewer, show data
%   E2_viewImage2D                      load 2D US data, initialize viewer and show data
%   E2_viewImage3D                      load 3D MRI data, initialize viewer and show data
%------------------------------------------------------------------------------
% IMAGE MODELS, INTERPOLATION, APPROXIMATION, SCALE, MULTILEVEL
%   E3_splineInterpolation2D            spline interpolation in 2D, from the book
%   E3_1D_basics                        basis concepts of interpolation in 1D
%   E3_1D_derivatives                   basis concepts of image models and its derivatives in 1D
%   E3_1D_scale                         concept of scale for 1D data
%   E3_2D_basics                        basis concepts of image models in 2D
%   E3_2D_derivative                    basis concepts of image models and its derivatives in 1D
%   E3_2D_generic                       use of interpolation and visualization,
%                                       loads data, initializes viewer and image model, shows results
%   E3_2D_scale                         concept of scale for 2D data
%   E3_getCellCenteredGrid              example for grid generation
%   E3_ij2xy3d                          example for 3D data indexing
%   E3_Hands_ij2xy                      coordinate transfer: (i,j) to (x,y))
%   E3_multilevel                       creates and visualizes multi-level representaion
%   E3_US_getMultilevel                 computes a multi-level representation of US images
%
%   E3_matlabInterpolation1D            matlab build in interpolator
%   E3_linearInterpolation1D            linear interpolation in 1D
%   E3_linearInterpolation2D            linear interpolation in 2D
%   E3_bsplines                         plots a B-spline
%   E3_checkDerivative                  performs a derivative check on splineInter
%   E3_splineInterpolation1D            spline interpolation in 1D, from the book
%   E3_truncatedSplineInterpolation1D   truncated spline interpolation in 1D
%   E3_MS_splineInterpolation1D         multi-scale spline model for 1D data
%   E3_MS_splineInterpolation2D         multi-scale spline model for 2D US data
%   E3_MS_splineInterpolation2Dext      multi-scale spline model for 2D US data  (more theta's)
%   E3_interpolation2D.m                interpolation in 2D (US data), various grid
%------------------------------------------------------------------------------
% IMAGE TRANSFORMATION
%   E4_Affine2D                         2D interpolation and affine transformations
%   E4_Affine2Dplain                    2D interpolation and affine transformations
%   E4_Bizarr                           2D interpolation and holomorphic transformations
%   E4_Rigid2D                          rotation of a 2D US image
%   E4_Rigid2Dplain                     rotation of a 2D US image, based on the book's version
%   E4_SplineTransformation2D           spline transformation of a 2D US image
%   E4_Translation2D                    translation of a 2D US image
%   E4_US_rotation                      rotating a 2D US image
%   E4_US_trafo                         2D interpolation and transformations of US image
%   E4_US_trafos                        various 2D interpolation and transformations of US image
%------------------------------------------------------------------------------
% LANDMARK BASED REGISTRATION
%   E5_2D_affine                        landmark based registration, affine transformation
%   E5_2D_quadratic                     landmark based registration, quadratic transformation
%   E5_Hands_TPS                        landmark based registration, thin-plate-spline transformation
%   E5_2D_TPS.m                         landmark based registration, tps transformation (various theta's)
%   E5_linear                           linear landmark based registration for hand example
%   E5_quadratic                        quadratic lm registration for hand example
%   E5_TPS                              thin plate spline lm registration for hand example
%   P5_LM                               generates nice plots
%------------------------------------------------------------------------------
% PIR - PARAMETRIC IMAGE REGISTRATION, MLPIR-MULTI_LEVEL PARAMETRIC IMAGE REGISTRATION
%   E6_quadrature_Spline1D              midpoint quadrature rule for 1D spline
%   E6_quadrature_Spline2D              midpoint quadrature rule for 2D spline
%   E6_quadrature_Gaussian2D            midpoint quadrature rule for 2D Gaussian
%   E6_quadrature_SSD2D                 SSD versus discretization width (midpoint quadrature rule)
%                                           
%   E6_FAIRplots                        demonstrates the use of FAIRplots.m
%                                           
%   E6_HNSP_SSD_rotation2D_level4       SSD versus rotations, linearInter, HNSP data, level=4
%   E6_HNSP_SSD_rotation2D_level8       SSD versus rotations, linearInter, HNSP data, level=8
%   E6_HNSP_SSD_translation2D_level4    SSD versus translations, linearInter, HNSP data, level=4
%   E6_HNSP_SSD_translation2D_level8    SSD versus translations, linearInter, HNSP data, level=4    
%   E6_HNSP_SSD_translation2D_level4_spline SSD versus translations, splineInter, HNSP data, level=4
%   E6_HNSP_SSD_translation2D_level8_spline SSD versus translations, splineInter, HNSP data, level=8
%   E6_HNSP_PIR_SSD_rotation2D_level4   SSD versus rotations, splineInter, HNSP data, level=4 
%   E6_HNSP_PIR_SSD_rotation2D_level7   SSD versus rotations, splineInter, HNSP data, level=7
%   E6_HNSP_PIR_SSD_rigid2D_level4      SSD versus rigid transformation, splineInter, HNSP data, level=4
%   E6_HNSP_PIR_SSD_rigid2D_level7      SSD versus rigid transformation, splineInter, HNSP data, level=7
%   E6_HNSP_PIR_SSD_affine2D_level5     SSD versus affine transformation, splineInter, HNSP data, level=5
%   E6_HNSP_PIR_SSD_spline2D_level5     SSD versus spline transformation, splineInter, HNSP data, level=5
%                                       
%   E6_Hands_PIR_Nesterov               PIR for 2D hand data with Nesterov optimization
%   E6_Hands_PIR_GN                     PIR for 2D hand data with Gauss-Newton optimization
%   E6_Hands_PIR_SD                     PIR for 2D hand data with steepest descent optimization
%   E6_Hands_affine                     ???
%   E6_HNSP_PIR_NelderMead              PIR for 2D HNSP data with Gauss-Newton optimization
%   E6_HNSP_PIR_GN                      PIR for 2D HNSP data with Gauss-Newton optimization
%                                       
%   E6_HNSP_RPIR                        regularized PIR, regularizer is a spline-moment matrix
%   E6_HNSP_PIR_scale                   MLPIR using a scale space
%
%   E6_Hands_MLPIR                      MLPIR for 2D hand data, straight
%   E6_PETCT_MLPIR                      MLPIR for 2D PET/CT data
%   E6_HNSP_MLPIR_SSD_rotation2D        MLPIR for 2D HNSP data with rotation, straight
%   E6_HNSP_MLPIR_SSD_rigid2D           MLPIR for 2D HNSP data  with rigid transformation, straight
%   E6_HNSP_MLPIR_SSD_affine2D          MLPIR for 2D HNSP data  with affine transformation, straight
%   E6_HNSP_MLPIR_reg                   MLPIR for 2D HNSP data  with splineTransformation2D, regularized
%   E6_3Dbrain_MLPIR_sparse             MLPIR for 3D brain data, uses splineTransformation3Dsparse
%------------------------------------------------------------------------------
% DISTANCE MEASURES
%   E7_histogram1D                      density estimation using a histogram
%   E7_histogram1D_ext                  density estimation using a histogram, extended version
%   E7_SSDforces                        distances and SSD forces for hand data
%   E7_HNSP_SSD_forces                  SSD and Force Fields
%                                       
%   E7_Hands_SSDvsRotation              various distances versus rotation angle
%   E7_PETCT_SSDvsRotation              various distances and Multi-Level Parametric Image Registration
%   E7_PETCT_MIvsRotation               mutual information versus rotation, 2D PET/CT data
%   E7_US_MIvsRotation                  mutual information versus rotation, 2D US data
%   E7_Hands_distance_rotation          various distances versus rotation angle, 2D hand data
%   E7_Hands_distance_rotation_ext      various distances versus rotation angle, 2D hand data, extended version
%   E7                                  various distances and Multi-Level Parametric Image Registration
%   E7_basic                            various distances measures for 2D PET/CT data
%   E7_extended                         various distances measures for 2D PET/CT data, extended version
%                                       
%   E7_PETCT_MLPIR                      various distances and MLPIR for 2D PET/CT data
%   E7_PETCT_MLPIR_ext                  various distances and MLPIR for 2D PET/CT data, extended version
%------------------------------------------------------------------------------
% REGULARIZATION
%   E8_elastic                          elastic, matrix based, matrix free
%   E8_curvature                        curvature, matrix based, matrix free
%   E8_forcesElastic                    deformation field and elastic forces
%   E8_forcesCurvature                  deformation field and curvature forces
%   E8_regularizationElasticMB          elastic regularization, matrix based
%   E8_regularizationElasticMF          elastic regularization, matrix free
%   E8_hyperElasticRegularizationMB     hyper elastic regularization, matrix based
%   E8_hyperElasticRegularizationMF     hyper elastic regularization, matrix based
%------------------------------------------------------------------------------
% NPIR - NON-PARAMETRIC IMAGE REGISTRATION and MLIR MULTILEVEL IMAGE REGISTRATION                                       
%   E9_HNSP_NPIR                        2D HNSP data, linearInter, SSD, mbElastic, GaussNewton
%   E9_HNSP_NPIR_pre                    2D HNSP data, linearInter, SSD, mbElastic, GaussNewton
%                                       with parametric pre-registration (rigid2D)
%   E9_Hands_NPIR                       2D hand data, splineInter, SSD, mbElastic, GaussNewton
%   E9_Hands_NPIR_MI_mbElas_BFGS        2D hand data, splineInter, MI, mbElastic, lBFGS
%   E9_Hands_NPIR_pre                   2D hand data, splineInter, SSD, mbElastic, GaussNewton
%                                       with parametric pre-registration (affine2D)
%   E9_Hands_NPIRmb_GN                  2D hand data, splineInter, MI, mbElastic, GaussNewton
%   E9_Hands_NPIRmf_GN                  2D hand data, splineInter, MI, mfElastic, GaussNewton
%   E9_Hands_NPIRmf_TR_nopre            2D hand data, splineInter, SSD, mfElastic, TrustRegion
%   E9_Hands_NPIRmf_TR_pcg              2D hand data, splineInter, SSD, mfElastic, TrustRegion
%
%   E9_Hands_NPIR_OPT                   checks a variety of optimizers
%     
%   E9_Hands_MLIR_SSD_mbElas            MLIR, 2D hand data, splineInter, SSD, mbElastic, GaussNewton
%   E9_Hands_MLIR_SSD_mfElas            MLIR, 2D hand data, splineInter, SSD, mfElastic, GaussNewton
%   E9_Hands_MLIR_SSD_mbCurv            MLIR, 2D hand data, splineInter, SSD, mbCurvature, GaussNewton
%   E9_Hands_MLIR_SSD_mfCurv            MLIR, 2D hand data, splineInter, SSD, mbCurvature, GaussNewton
%   E9_Hands_MLIR_wSSD_mfElas           MLIR, 2D hand data, splineInter, weighted SSD, mbCurvature, GaussNewton
%   E9_HNSP_MLIR_SSD_mbElas             MLIR, 2D HNSP data, splineInter, SSD, mbElastic, GaussNewton
%   E9_HNSP_MLIR_SSD_mfElas             MLIR, 2D HNSP data, splineInter, SSD, mfElastic, GaussNewton
%   E9_HNSP_MLIR_SSD_mbCurv             MLIR, 2D HNSP data, splineInter, SSD, mbCurvature, GaussNewton
%   E9_HNSP_MLIR_SSD_mfCurv             MLIR, 2D HNSP data, splineInter, SSD, mbCurvature, GaussNewton
%    
%   E9_MRIhead_MLIR_SSD_mbElas          MLIR, 2D MRI head data, splineInter, SSD, mbElastic, GaussNewton
%   E9_MRIhead_MLIR_MI_mbElas           MLIR, 2D MRI head data, splineInter, MI,  mbElastic, GaussNewton
%   E9_MRIhead_MLIR_NGF_mbElas          MLIR, 2D MRI head data, splineInter, NGF, mbElastic, GaussNewton
%   E9_MRIhead_MLIRlBFGS_MI_mfElas      MLIR, 2D MRI head data, splineInter, MI,  mbElastic, lBFGS
%   E9_MRIhead_MLIRlBFGS_NGF_mbElas     MLIR, 2D MRI head data, splineInter, NGF, mbElastic, lBFGS
%   E9_PETCT_MLIR_NGF_mbElas            MLIR, 2D PET/CT data,   splineInter, NGF, mbElastic, GaussNewton
%   E9_PETCT_MLIR_NGF_mfElas            MLIR, 2D PET/CT data,   splineInter, NGF, mfElastic, GaussNewton
%   E9_PETCT_MLIR_NGF_mbCurv            MLIR, 2D PET/CT data,   splineInter, NGF, mbCurvature, GaussNewton
%   E9_PETCT_MLIR_NGF_mfCurv            MLIR, 2D PET/CT data,   splineInter, NGF, mfCurvature, GaussNewton
%   E9_PETCT_MLIRlBFGS_NGF_mbElas       MLIR, 2D PET/CT data,   splineInter, NGF, mbElastic, lBFGS
%   E9_PETCT_MLIRlBFGS_NGF_mbCurv       MLIR, 2D PET/CT data,   splineInter, NGF, mbCurvature, lBFGS
%   E9_PETCT_MLIRlBFGS_MI_mbElas        MLIR, 2D PET/CT data,   splineInter, MI,  mbElastic, lBFGS
%   E9_PETCT_MLIRlBFGS_MI_mfElas        MLIR, 2D PET/CT data,   splineInter, MI,  mfElastic, lBFGS
%   E9_PETCT_MLIRlBFGS_MI_mbCurv        MLIR, 2D PET/CT data,   splineInter, MI,  mbCurvature, lBFGS
%   E9_PETCT_MLIRlBFGS_MI_mfCurv        MLIR, 2D PET/CT data,   splineInter, MI,  mfCurvature, lBFGS
%                  
%   E9_Hands_MSIR                       Multi-scale registration, 
%                                       2D hand data, splineInter, SSD, mbCurvature, GaussNewton
%   E9_HNSP_MLIR_TR                     MLIR, 2D hand data, splineInter, SSD, mbElastic, Trust Region
%
%   E9_3Dbrain_GN                       MLIR, 3D brain data, splineInter, SSD, mfElastic, GaussNewton
%   E9_3Dbrain_TR                       MLIR, 3D brain data, splineInter, SSD, mfElastic, TrustRegion
%   E9_3Dbrain_lBFGS                    MLIR, 3D brain data, splineInter, SSD, mfElastic, lBFGS
%         
% HYPERELASTICITY
%   Ehyper_2Ddisc2C_MB                  Hyper elasticity, 2D, matrix based          
%   Ehyper_3Dbrain_MF                   Hyper elasticity, 3D, matrix free          
%==============================================================================

function debit = contents
if nargout == 0, 
    help(mfilename); 
	return; 
end;

debit = {
    'contents.m'
    'spline1D.m'
    'testExamples.m'


    'E2_setup2DhandData.m'
    'E2_setup3DbrainData.m'
    'E2_viewImage2D.m'
    'E2_viewImage3D.m'
    

    'E3_splineInterpolation2D.m'
    'E3_1D_basics.m'
    'E3_1D_derivatives.m'
    'E3_1D_scale.m'
    'E3_2D_basics.m'
    'E3_2D_derivative.m'
    'E3_2D_generic.m'
    'E3_2D_scale.m'
    'E3_getCellCenteredGrid.m'
    'E3_ij2xy3d.m'
    'E3_Hands_ij2xy.m'
    'E3_multilevel.m'
    'E3_US_getMultilevel.m'
    'E3_matlabInterpolation1D.m'
    'E3_linearInterpolation1D.m'
    'E3_linearInterpolation2D.m'
    'E3_bsplines.m'
    'E3_checkDerivative.m'
    'E3_splineInterpolation1D.m'
    'E3_truncatedSplineInterpolation1D.m'
    'E3_MS_splineInterpolation1D.m'
    'E3_MS_splineInterpolation2D.m'
    'E3_MS_splineInterpolation2Dext.m'
    'E3_interpolation2D.m'
        

    'E4_Affine2D.m'
    'E4_Affine2Dplain.m'
    'E4_Bizarr.m'
    'E4_Rigid2D.m'
    'E4_Rigid2Dplain.m'
    'E4_SplineTransformation2D.m'
    'E4_Translation2D.m'
    'E4_US_rotation.m' 
    'E4_US_trafo.m'
    'E4_US_trafos.m'
    

    'E5_2D_affine.m'
    'E5_2D_quadratic.m'
    'E5_Hands_TPS.m'
    'E5_2D_TPS.m'  
    'E5_linear.m'
    'E5_quadratic.m'
    'E5_TPS.m'
    'P5_LM.m'
    

    'E6_quadrature_Spline1D.m'
    'E6_quadrature_Spline2D.m'
    'E6_quadrature_Gaussian2D.m'
    'E6_quadrature_SSD2D.m'
    
    'E6_FAIRplots.m'

	'E6_HNSP_SSD_rotation2D_level4.m'
    'E6_HNSP_SSD_rotation2D_level8.m'
    'E6_HNSP_SSD_translation2D_level4.m'
    'E6_HNSP_SSD_translation2D_level8.m'
    'E6_HNSP_SSD_translation2D_level4_spline.m'
    'E6_HNSP_SSD_translation2D_level8_spline.m'
    'E6_HNSP_PIR_SSD_rotation2D_level4.m'
    'E6_HNSP_PIR_SSD_rotation2D_level7.m'
    'E6_HNSP_PIR_SSD_rigid2D_level4.m'
    'E6_HNSP_PIR_SSD_rigid2D_level7.m'
    'E6_HNSP_PIR_SSD_affine2D_level5.m'
    'E6_HNSP_PIR_SSD_spline2D_level5.m'

    'E6_Hands_PIR_Nesterov.m'
    'E6_Hands_PIR_GN.m'
    'E6_Hands_PIR_SD.m'
    'E6_Hands_affine.m'
    'E6_HNSP_PIR_NelderMead.m'
    'E6_HNSP_PIR_GN.m'
    'E6_HNSP_RPIR.m'
    'E6_HNSP_PIR_scale.m'

    'E6_Hands_MLPIR.m'
    'E6_PETCT_MLPIR.m'
    'E6_HNSP_MLPIR_SSD_rotation2D.m'
    'E6_HNSP_MLPIR_SSD_rigid2D.m'
    'E6_HNSP_MLPIR_SSD_affine2D.m'
    'E6_HNSP_MLPIR_reg.m'
    'E6_3Dbrain_MLPIR_sparse.m'


    'E7_histogram1D.m'
    'E7_histogram1D_ext.m'
    'E7_SSDforces.m'
    'E7_HNSP_SSD_forces.m'

    'E7_Hands_SSDvsRotation.m'
    'E7_PETCT_SSDvsRotation.m'
    'E7_PETCT_MIvsRotation.m'
    'E7_US_MIvsRotation.m'
    'E7_Hands_distance_rotation.m'
    'E7_Hands_distance_rotation_ext.m'
    'E7.m'
    'E7_basic.m'
    'E7_extended.m'

    'E7_PETCT_MLPIR.m'
    'E7_PETCT_MLPIR_ext.m'


    'E8_elastic.m'
    'E8_curvature.m'
    'E8_forcesElastic.m'
    'E8_forcesCurvature.m'
    'E8_regularizationElasticMB.m'
    'E8_regularizationElasticMF.m'
    'E8_RegularizationHyperElasticMB.m'
    'E8_RegularizationHyperElasticMF.m'


    'E9_HNSP_NPIR.m'
    'E9_HNSP_NPIR_pre.m'
    'E9_Hands_NPIR.m'
    'E9_Hands_NPIR_MI_mbElas_BFGS.m'
    'E9_Hands_NPIR_pre.m'
    'E9_Hands_NPIRmb_GN.m'
    'E9_Hands_NPIRmf_GN.m'
    'E9_Hands_NPIRmf_TR_nopre.m'  
    'E9_Hands_NPIRmf_TR_pcg.m'   
	
	'E9_Hands_NPIR_OPT.m'
   
    'E9_Hands_MLIR_SSD_mbElas.m'
    'E9_Hands_MLIR_SSD_mfElas.m'
    'E9_Hands_MLIR_SSD_mbCurv.m'
    'E9_Hands_MLIR_SSD_mfCurv.m'
    'E9_Hands_MLIR_wSSD_mfElas.m'
    'E9_HNSP_MLIR_SSD_mbElas.m'
    'E9_HNSP_MLIR_SSD_mfElas.m'
    'E9_HNSP_MLIR_SSD_mbCurv.m'
    'E9_HNSP_MLIR_SSD_mfCurv.m'
    
    'E9_MRIhead_MLIR_SSD_mbElas.m'
    'E9_MRIhead_MLIR_MI_mbElas.m'
    'E9_MRIhead_MLIR_NGF_mbElas.m'
    'E9_MRIhead_MLIRlBFGS_MI_mfElas.m'
    'E9_MRIhead_MLIRlBFGS_NGF_mbElas.m'
	
    'E9_PETCT_MLIR_NGF_mbElas.m'
    'E9_PETCT_MLIR_NGF_mfElas.m'
    'E9_PETCT_MLIR_NGF_mbCurv.m'
    'E9_PETCT_MLIR_NGF_mfCurv.m'
	  'E9_PETCT_MLIRlBFGS_NGF_mbElas.m'
    'E9_PETCT_MLIRlBFGS_NGF_mbCurv.m'
	  'E9_PETCT_MLIRlBFGS_MI_mbElas.m'
    'E9_PETCT_MLIRlBFGS_MI_mfElas.m'
    'E9_PETCT_MLIRlBFGS_MI_mbCurv.m'
    'E9_PETCT_MLIRlBFGS_MI_mfCurv.m'

    'E9_Hands_MSIR.m'
    'E9_HNSP_MLIR_TR.m'
	
    'E9_3Dbrain_GN.m'
    'E9_3Dbrain_TR.m'
    'E9_3Dbrain_lBFGS.m'


    'Ehyper_2Ddisc2C_MB.m'
	'Ehyper_3Dbrain_MF.m'		
    };
%==============================================================================
