% ==================================================================================
% (c) Lars Ruthotto 2012/07/26, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/lars-ruthotto.html
%
% Mass-Preserving Registration of 3D PET of a mice heart using FEM
%
% ==================================================================================
close all; clc; clear;
%setup3DmiceData; 
setup3DmouseData;

level = 3;
omega = ML{level}.omega; m = ML{level}.m;

viewImage('reset','viewImage','imgmontage','colormap','gray');
imgModel('reset','imgModel','splineInterMex','regularizer','moments','theta',1e-1);
regularizer('reset','regularizer','mbHyperElasticFEM','alpha',1e1);
[T,R] = imgModel('coefficients',ML{level}.T,ML{level}.R,omega,'out',0);

Mesh = TetraMesh1(omega,m);
xc   = Mesh.xn;
Rc = imgModel(R,omega,Mesh.mfPi(xc,'C'));
Tc = imgModel(T,omega,Mesh.mfPi(xc,'C'));
fctn = @(yc) FEMPIREobjFctn(T,Rc,Mesh,xc,yc(:));
fctn([]);

FAIRplotsFEM('reset','mode','FEM','fig',level,'plots',1);
FAIRplotsFEM('init',struct('Tc',T,'Rc',R,'Mesh',Mesh));

yc = GaussNewton(fctn,xc(:),'Plots',@FAIRplotsFEM,'solver','mbPCG-Jacobi');

% matrixfree version
regularizer('reset','regularizer','mfHyperElasticFEM','alpha',1e1);
ycMF = GaussNewton(fctn,xc(:),'Plots',@FAIRplotsFEM,'solver',@FEMPIREsolveGN_PCG);

%showResults(ML,yc,'level',level);
