% ==================================================================================
% (c) Lars Ruthotto 2011/02/08, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/lars-ruthotto.html
%
% 2D example for hyperelastic image registration using a finite element
% method and a uniform refinement strategy
%
% ==================================================================================

setup2DhandData
trafo('reset','trafo','rigid2D');
distance('reset','distance','SSD');
imgModel('reset','imgModel','splineInter','regularizer','moments','theta',4e-1);
regularizer('reset','regularizer','mbElasticFEM','alpha',1e2);
maxLevel = 6;

[yOpt,wc,his] = MLIRFEM(ML,'minLevel',4,'maxLevel',maxLevel,'parametric',0);
