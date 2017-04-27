% ==================================================================================
% (c) Lars Ruthotto 2012/07/26, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/lars-ruthotto.html
%
% 2D Mass-Preserving Registration of Gaussian Blobs
%
% ==================================================================================
close all; clc; clear;
setup2DGaussianData;

distance('reset','distance','SSD');
imgModel('reset','imgModel','splineInter','regularizer','moments','theta',4e-1);
regularizer('reset','regularizer','mbHyperElasticFEM','alpha',1e2,'alphaLength',1,'alphaVolume',1);

NPIRpara    = optPara('NPIR-GN');
NPIRpara.Plots = @FAIRplotsFEM;
NPIRpara.lineSearch = @ArmijoDiffeomorphicFEM;

yc = MLIRFEM(ML,'parametric',0,'NPIRpara',NPIRpara,'NPIRobj',@FEMPIREobjFctn);
