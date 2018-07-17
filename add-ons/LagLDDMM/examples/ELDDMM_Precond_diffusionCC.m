% =========================================================================
% This code is part of the Matlab-based toolbox 
% LagLDDDM - A Lagrangian Gauss--Newton--Krylov Solver for Mass- and 
%                        Intensity-Preserving Diffeomorphic Image Registration
% 
% For details and license info see 
% - https://github.com/C4IR/FAIR.m/tree/master/add-ons/LagLDDMM
%
% Example illustrating the performance of different preconditioners for
% the Gauss-Newton system in LDDMM. The example is described in detail in
% Section 4.2 of the paper:
%
% @article{MangRuthotto2017,
%   Title = {A {L}agrangian {G}auss--{N}ewton--{K}rylov solver for mass- and intensity-preserving diffeomorphic image registration},
%   Year = {2017},
%   Journal = {SIAM Journal on Scientific Computing},
%   Author = {A. Mang, L. Ruthotto},
% }
%
% =========================================================================
clear all; clc;
setup2Ddisc2CData

%% run affine pre-registration
imgModel('reset','imgModel','splineInterMex','regularizer','moments','theta',.1);
        trafo('set','trafo','affine2D');
        distance('set','distance','SSD');

alpha = [6e2 0];
parametric = 0;
pad  = 0.5;
nt = 0;
N    = 3;
mV     = @(m) ceil(1*m);
minLevel = 5;
maxLevel = 5;

%% run multilevel LDDMM registration

% 1) setup grid for velocities (padded)
omegaV = omega; omegaV(1:2:end) = omegaV(1:2:end)-pad;  omegaV(2:2:end) = omega(2:2:end)+pad;

% 2) setup regularizer
regularizer('reset','regularizer','mfDiffusionCC','alpha',alpha,'nt',nt,'HessianShift',1e-2); % stationary velocity

NPIRpara    = optPara('NPIR-GN');
NPIRpara.maxIter = 40;
NPIRpara.scheme = @GaussNewtonLDDMM;
[vc,~,wc,his] = MLLDDMM(ML,'minLevel',minLevel,'maxLevel',maxLevel,...
    'omegaV',omegaV,'mV',mV,'N',N,'parametric',parametric,'NPIRpara',NPIRpara,'plots',1);


%% generate plots to explore conditioning / sparsity
xc = getNodalGrid(omega,ML{minLevel}.m);
[T,R] = imgModel('coefficients',ML{minLevel}.T,ML{minLevel}.R,omega);
Rc    = imgModel(R,omega,center(xc,ML{minLevel}.m));
  
mV = ML{minLevel}.m;
fctn = @(vc) LDDMMobjFctn(T,Rc,omega,ML{minLevel}.m,0*vc,center(xc,ML{minLevel}.m),omegaV,mV,N,vc);

close all
regularizer('set','regularizer','mbDiffusionCC')
[Sc,dS,d2Smb] = regularizer(vc,omegaV,mV);
% get objective function
[J0,p0,dJ0,Hmb0] = fctn(0*vc);
[Jc,pc,dJc,Hmbc] = fctn(vc);

regularizer('set','regularizer','mfDiffusionCC')
[J0,p0,dJ0,Hmf0] = fctn(0*vc);
[Jc,pc,dJc,Hmfc] = fctn(vc);

%%
fig = figure(); clf;
fig.Name =sprintf('%s: Sparsity of Hessian',mfilename);

subplot(1,2,1)
spy(Hmb0);
title('H(v0)');

subplot(1,2,2);
spy(Hmbc);
title('H(vOpt');
%% explore PCG convergence at first iteration
hd = prod((omega(2:2:end)-omega(1:2:end))./ML{minLevel}.m);
Hmb0 = Hmb0 + hd* regularizer('get','HessianShift')* speye(size(Hmb0));
D   = diag(Hmb0); % D is diagonal
L   = tril(Hmb0); % Symmetric Gauss Seidel Preconditioning,

PCjac = @(x) D.\x;
PCsgs = @(x) full(L\(D.*(L'\x)));
PCnone = @(x) x;

[~,iter0cg,relres0cg,resvec0cg] = spectralPrecondPCG(-dJ0(:), Hmf0, 250, 1e-10,'prec',PCnone);
[~,iter0jac,relres0jac,resvec0jac] = spectralPrecondPCG(-dJ0(:), Hmf0, 250, 1e-10,'prec',PCjac);
[~,iter0sgs,relres0sgs,resvec0sgs] = spectralPrecondPCG(-dJ0(:), Hmf0, 250, 1e-10,'prec',PCsgs);
[~,iter0sp,relres0sp,resvec0sp] = spectralPrecondPCG(-dJ0(:), Hmf0, 250, 1e-10);

%%
figConv = figure(); clf
figConv.Name = 'PCG convergence';

subplot(1,2,1);
semilogy(resvec0cg/resvec0cg(1));
hold on;
semilogy(resvec0jac/resvec0jac(1));
semilogy(resvec0sgs/resvec0sgs(1));
semilogy(resvec0sp/resvec0sp(1));
title('first GN iteration');
ax = [0 250 1e-10 10];
axis(ax)
legend('CG','PCG-Jac','PCG-SGS','PCG-Spec','Location','SouthWest')


%% explore PCG convergence at final iteration
hd = prod((omega(2:2:end)-omega(1:2:end))./ML{minLevel}.m);
Hmb0 = Hmb0 + hd* regularizer('get','HessianShift')* speye(size(Hmb0));
D   = diag(Hmbc); % D is diagonal
L   = tril(Hmbc); % Symmetric Gauss Seidel Preconditioning,

linSol = @spectralPrecondPCG;
PCjac = @(x) D.\x;
PCsgs = @(x) full(L\(D.*(L'\x)));
PCnone = @(x) x;

[~,itercg,relrescg,resveccg]    = linSol(-dJc(:), Hmfc, 250, 1e-10,'prec',PCnone);
[~,iterjac,relresjac,resvecjac] = linSol(-dJc(:), Hmfc, 250, 1e-10,'prec',PCjac);
[~,itersgs,relressgs,resvecsgs] = linSol(-dJc(:), Hmfc, 250, 1e-10,'prec',PCsgs);
[~,itersp,relressp,resvecsp]    = linSol(-dJc(:), Hmfc, 250, 1e-10);

%%
fig = figure(figConv); 
subplot(1,2,2);
semilogy(resveccg/resveccg(1));
hold on;
semilogy(resvecjac/resvecjac(1));
semilogy(resvecsgs/resvecsgs(1));
semilogy(resvecsp/resvecsp(1));
title('final iteration');

legend('CG','PCG-Jac','PCG-SGS','PCG-Spec','Location','SouthWest')
axis(ax)

