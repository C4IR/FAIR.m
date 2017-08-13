% ==================================================================================
% (c) Lars Ruthotto 2012/07/26, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/lars-ruthotto.html
%
% Verification code for matrixfree implementation
%
% 'relerr' shows the relative difference between matrixfree and matrix
% based implementations
%
% ==================================================================================
close all; clc; clear;
setup2DGaussianData;

distance('reset','distance','SSD');
imgModel('reset','imgModel','splineInterMex','regularizer','moments','theta',4e-1);
regularizer('reset','regularizer','mbHyperElasticFEM','alpha',5e1);

minLevel = 3;
maxLevel = 8;

NPIRpara            = optPara('NPIR-GN');
NPIRpara.Plots      = @FAIRplotsFEM;
NPIRpara.lineSearch = @ArmijoDiffeomorphicFEM;
NPIRpara.solver     = 'mbPCG-Jacobi';

[yMP_MB,~,hisMP_MB] = MLIRFEM(ML,'parametric',0,...
      'NPIRobj',@FEMPIREobjFctn,'NPIRpara',NPIRpara,...
      'minLevel',minLevel,'maxLevel',maxLevel);

regularizer('reset','regularizer','mfHyperElasticFEM','alpha',5e1);
NPIRpara.solver = @FEMPIREsolveGN_PCG;%'PCG-hyperElastic';

[yMP_MF,~,hisMP_MF] = MLIRFEM(ML,'parametric',0,...
    'NPIRobj',@FEMPIREobjFctn,'NPIRpara',NPIRpara,...
    'minLevel',minLevel,'maxLevel',maxLevel);

relerr = norm(yMP_MB-yMP_MF)/norm(yMP_MB)

%%
%save('EFEM_SSDvsMP','ML','yFEMPIRE','ySSD','hisFEMPIRE','hisSSD');