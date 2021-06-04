%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
% This is a testing environment for the files in the folder
% kernel/regularizers
% 1. Based on regularizers/contents.m, a list of required files is generated and it 
%    is verified, that all files are present; additional files are listed.
% 2. All c-files are compiled.
% 3. All files are executed.
% 4. Check administrative module regularizer.m
% 5. Test specific regularizer 
%==============================================================================

FAIRcheckFiles(mfilename);

%% 4. Test module regularizer,
% setup dummy parameters
alpha  = 1;

% for (linear) elastic regularizer
mu     = 1;
lambda = 0;

% for hyperelastic regularizer
alphaLength = 13.2;
alphaArea   =  1.3;
alphaVolume =  3.4;

% linear elastic
regularizer('reset','regularizer','mfElastic',...
  'alpha',alpha,'mu',mu,'lambda',lambda);
regularizer('disp');
[scheme,parameter] = regularizer;
parameter = cell2struct({parameter{2:2:end}},{parameter{1:2:end}},2)
regularizer('clear');

% curvature
regularizer('reset','regularizer','mbCurvature','alpha',alpha);
regularizer('disp');
[scheme,parameter] = regularizer;
parameter = cell2struct({parameter{2:2:end}},{parameter{1:2:end}},2)
regularizer('clear');

% hyperelastic
regularizer('reset','regularizer','mfHyperElastic','alpha',alpha,...
    'alphaLength',alphaLength,'alphaArea',alphaArea,'alphaVolume',alphaVolume);
regularizer('disp');
[scheme,parameter] = regularizer;
parameter = cell2struct({parameter{2:2:end}},{parameter{1:2:end}},2)
regularizer('clear');

%% Test mb versus mf and derivatives
regularizer('reset','regularizer','mfElastic',...
  'alpha',alpha,'mu',mu,'lambda',lambda);
regularizer('disp')

omega = [0,1,0,1];
m     = [32,32]/8;
hd    = prod((omega(2:2:end)-omega(1:2:end))./m);

% check elasticity
H.omega  = omega; 
H.m      = m; 
H.alpha  = regularizer('get','alpha');
H.mu     = regularizer('get','mu');
H.lambda = regularizer('get','lambda');
H.regularizer = regularizer;

X   = getStaggeredGrid(omega,m);
B  = getElasticMatrixStg(omega,m,mu,lambda);
Y  = randn(size(X));

M   = speye(length(Y),length(Y));
A   = M + hd* H.alpha *B'*B;
rhs = A*X;
% rhs = A(:,5);
u0  = zeros(size(Y));

% prepare for multigrid
H.MGlevel      = log2(m(1))+1;
H.MGcycle      = 4;
H.MGomega      = 2/3;   %% !!! 0.5 should be better
H.MGsmoother   = 'mfJacobi';
H.MGpresmooth  = 10;
H.MGpostsmooth = 10;
H.d2D.M        = full(diag(M));

[Sc, dS, d2S] = regularizer(X,omega,m);
H.d2S = d2S;

Zmg = mfvcycle(H,u0,rhs,1e-12,H.MGlevel,5);
testMG = norm( A\rhs-Zmg);
fprintf('|A\\rhs-uMG|=%s\n',num2str(testMG));
assert(testMG<1e-10,...
  'significant difference between matrix-based and matrix free code');

%% check matrix free implementation
Omega = {0,[0,1,0,2],[0,1,0,2,0,3]}
M     = {0,[4,5],[3,4,5]}

%%

Regs = {'Elastic','Curvature','HyperElastic'};

Para = {
  {'alpha',alpha,'mu',mu,'lambda',lambda}
  {'alpha',alpha}
  {,'alpha',alpha,'alphaLength',alphaLength,'alphaArea',alphaArea,...
    'alphaVolume',alphaVolume}
};

Grids = {
  @(omega,m) getStaggeredGrid(omega,m)
  @(omega,m) getCellCenteredGrid(omega,m)
  @(omega,m) getNodalGrid(omega,m)
};

relErr = @(a,b) norm(a - b) / (norm(b) + (norm(b)==0));

for k=1:length(Regs),
  for dim = 2:3,
    
    % test mb versus mf
    fprintf('test mb%s versus mf%s, spatial dimension is %d\n',...
      Regs{k},Regs{k},dim);
    
    omega = Omega{dim};
    m     = M{dim}
    X     = Grids{k}(omega,m);
    Y     = randn(size(X));
 
    regularizer('reset','regularizer',['mb',Regs{k}],Para{k}{:});
%     regularizer('disp')
    [mbS,mbdS,mbd2S] = regularizer(Y,omega,m);
    
    regularizer('reset','regularizer',['mf',Regs{k}],Para{k}{:});
%     regularizer('disp')
    [mfS,mfdS,mfd2S] = regularizer(Y,omega,m);
    
    e1 = relErr(mbS,mfS);
    fprintf('%-30s = %s\n','norm(mfS  - mbS )/norm(mbS)',  num2str(e1))
    e2 = relErr(mbdS,mfdS);
    fprintf('%-30s = %s\n','norm(mfdS - mbdS)/norm(mbdS)', num2str(e2))
 
    assert(e1+e2<1e-10,...
      'significant difference between matrix-based and matrix free code');
    
    % test dericatives
    fprintf('test derivative of [%s]\n',regularizer);
    
    fctn = @(Y) regularizer(Y,omega,m);
    checkDerivative(fctn,Y);
    title(sprintf('checkDerivative: %s - %dD',regularizer,dim));

  end;
end

testEnd;
%==============================================================================

