folder = fullfile(fileparts(which(mfilename)),'examples');
files = {
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
        };

FAIRcheckFiles(folder,files,mfilename,'FAIRkeyboard','off');    

 