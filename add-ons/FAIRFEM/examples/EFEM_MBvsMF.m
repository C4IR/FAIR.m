% ==================================================================================
% (c) Lars Ruthotto 2012/07/26, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/lars-ruthotto.html
%
% Simplified Registration of Cardiac PET
%
% data Gaussian Blobs
%
% ==================================================================================
close all; clc; clear;
setup2DGaussianData;

distance('reset','distance','SSD');
imgModel('reset','imgModel','splineInterMex','regularizer','moments','theta',4e-1);
regularizer('reset','regularizer','mbHyperElasticFEM','alpha',1e2);

minLevel = 3;
maxLevel = 3;

NPIRpara            = optPara('NPIR-GN');
NPIRpara.Plots      = @FAIRplotsFEM;
NPIRpara.lineSearch = @ArmijoDiffeomorphicFEM;
NPIRpara.solver     = 'mbPCG-Jacobi';

[ySSD_MB,~,hisSSD_MB] = MLIRFEM(ML,'parametric',0,...
    'NPIRobj',@FEMobjFctn,'NPIRpara',NPIRpara,...
    'minLevel',minLevel,'maxLevel',maxLevel);

regularizer('reset','regularizer','mfHyperElasticFEM','alpha',5e1);
NPIRpara.solver = @FEMSSDsolveGN_PCG;

[ySSD_MF,~,hisSSD_MF] = MLIRFEM(ML,'parametric',0,...
    'NPIRobj',@FEMobjFctn,'NPIRpara',NPIRpara,...
    'minLevel',minLevel,'maxLevel',maxLevel);

relerr = norm(ySSD_MB-ySSD_MF)/norm(ySSD_MB)

%%
%save('EFEM_SSDvsMP','ML','yFEMPIRE','ySSD','hisFEMPIRE','hisSSD');