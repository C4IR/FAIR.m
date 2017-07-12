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
regularizer('reset','regularizer','mbHyperElasticFEM','alpha',5e1);

minLevel = 5;
maxLevel = 7;

NPIRpara            = optPara('NPIR-GN');
NPIRpara.Plots      = @FAIRplotsFEM;
NPIRpara.lineSearch = @ArmijoDiffeomorphicFEM;
NPIRpara.solver     = 'mbPCG-Jacobi';


[yFEMPIRE,~,hisFEMPIRE] = MLIRFEM(ML,'parametric',0,...
    'NPIRobj',@FEMPIREobjFctn,'NPIRpara',NPIRpara,...
    'minLevel',minLevel,'maxLevel',maxLevel);

[ySSD,~,hisSSD] = MLIRFEM(ML,'parametric',0,...
    'NPIRobj',@FEMobjFctn,'NPIRpara',NPIRpara,...
    'minLevel',minLevel,'maxLevel',maxLevel);
%%
save('EFEM_SSDvsMP','ML','yFEMPIRE','ySSD','hisFEMPIRE','hisSSD');