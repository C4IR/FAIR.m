%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% Tutorial for FAIR: curvature regularizer, matrix based and matrix free
%
% illustrates the tools for L2-norm based regularization
%
%  S(y) = alpha/2 * norm(B*(y-yRef)^2,
%
% where
%  alpha regularization parameter, weights regularization versus 
%        distance in the joint objective function, alpha = 1 here
%  yRef  is a reference configuration, e.g. yRef = x or a 
%        pre-registration result
%  B     a discretized partial differential operator either in explicit
%        matrix form or as a structure containing the necessary
%        parameters to compute B*y, here B = curvature
%==============================================================================

clear, close all, help(mfilename);

%% 2D example
omega  = [0,1,0,1]; % physical domin
m      = [16,12];   % number of discretization points
hd     = prod((omega(2:2:end)-omega(1:2:end))./m);

alpha  = 1;         % regularization parameter, irrelevant in this tutorial

yRef   = getCellCenteredGrid(omega,m); % reference for regularization
yc     = rand(size(yRef));             % random transformation

%% 1. build curvature operator on a cell-centered grid 
%  - explicitly 
%  - via regularizer.m
%  - via A(:,j) = B(e_j), B is the matrix free version

B = getCurvatureMatrix(omega,m);
Amb = hd*alpha*(B'*B);

regularizer('reset','regularizer','mbCurvature','alpha',alpha);
% regularizer('disp');
[Smb,dSmb,d2Smb] = regularizer(yc-yRef,omega,m);

regularizer('reset','regularizer','mfCurvature','alpha',alpha);
% regularizer('disp');
[Scmf, dSmf, d2Smf] = regularizer(yc,omega,m); 

% note: B*y can also be computed without storing B. The following loop 
% computes A(:,j) = mfBy(e(j)), where mfBy is the matrix free operation B*y
% and e(j) is the j-th unit vector; hence A == B, 

e = @(j) ((1:size(B,2))' == j); % j-th unit vector
A = sparse(size(B,1),size(B,2));
for j=1:size(B,2), A(:,j) = d2Smf.By(e(j),omega,m); end;

% now the transpose
e = @(j) ((1:size(B,1))' == j); % j-th unit vector
C = sparse(size(B,2),size(B,1));
for j=1:size(B,1), C(:,j) = d2Smf.BTy(e(j),omega,m); end;

% visualize matrices and transpose
FAIRfigure(1); clf; 
subplot(2,3,1); spy(B);    title('B =curvature operator on cell-centered grid')
subplot(2,3,2); spy(A);    title('B from mfBy');
subplot(2,3,3); spy(B-A);  title('difference |B-mfBy()|');
subplot(2,3,4); spy(B');   title('B''')
subplot(2,3,5); spy(C);    title('B'' from mfBy');
subplot(2,3,6); spy(B'-C); title('difference |B''-mfBy(...,''BTy'')|');


%% 3D example
omega  = [0,4,0,2,0,1]; % physical domin
m      = [16,12,8];     % number of discretization points
hd     = prod((omega(2:2:end)-omega(1:2:end))./m);

yRef   = getCellCenteredGrid(omega,m); % reference for regularization
yc     = rand(size(yRef));             % random transformation

% build elasticity operator on a staggered grid
B = getCurvatureMatrix(omega,m);
Amb = hd*alpha*(B'*B);

regularizer('reset','regularizer','mbCurvature','alpha',alpha);
% regularizer('disp');
[Smb,dSmb,d2Smb] = regularizer(yc-yRef,omega,m);

regularizer('reset','regularizer','mfCurvature','alpha',alpha);
% regularizer('disp');
[Smf, dSmf, d2Smf] = regularizer(yc-yRef,omega,m); 


zc     = randn(size(B,1),1);        % random for the adjoint

FAIRfigure(2); clf; 
subplot(2,2,1); spy(B);          title('B elastic on staggered grid')
subplot(2,2,3); spy(Amb);        title('mb: hd\alpha B''*B')
subplot(2,2,4); spy(d2Smb);      title('mb: d2S')
subplot(2,2,2); spy(Amb-d2Smb);  title('difference')

Splain = 0.5*(yc-yRef)'*Amb*(yc-yRef);

fprintf('Splain = %-10s\n',num2str(Splain))
fprintf('Smb    = %-10s, Splain-Smb =%-10s\n',num2str(Smb),num2str(Splain-Smb))
fprintf('Smf    = %-10s, Splain-Smf =%-10s\n',num2str(Smf),num2str(Splain-Smf))

compare = @(s,l,r) fprintf('%-25s = %s\n',s,num2str(norm(l-r)));
compare('||B*yc-mfBy(y)||',B*yc,d2Smf.By(yc,omega,m));
compare('||B''*z-mfBy(z,''BTy'')||',B'*zc,d2Smf.BTy(zc,omega,m));


%% check curvature

%% 3D example
omega  = [0,4,0,2,0,1]; % physical domin
m      = [32,16,8];     % number of discretization points
hd     = prod((omega(2:2:end)-omega(1:2:end))./m);

yRef   = getCellCenteredGrid(omega,m); % reference for regularization
yc     = rand(size(yRef));             % random transformation

% build elasticity operator on a staggered grid
B = getCurvatureMatrix(omega,m);

regularizer('reset','regularizer','mfCurvature','alpha',alpha);
% regularizer('disp');
[Smf, dSmf, d2Smf] = regularizer(yc-yRef,omega,m); 

H.d2S = d2Smf;
H.omega  = omega;
H.m      = m;
H.alpha  = regularizer('get','alpha');
H.regularizer = regularizer;

By = B*yc;
testMF  = norm(By-d2Smf.By(yc,omega,m));
testMFT = norm(B'*By-d2Smf.BTy(By,omega,m));
fprintf('|By-mfBy|     = %s\n',num2str(testMF));
fprintf('|B''*By-mfBy''| = %s\n',num2str(testMFT));
% ---------------------------------------------------------------------
% using multigrid to solve A*uc = fc, where A = I + alpha*hd*B'*B

xc  = getCellCenteredGrid(omega,m);
M   = speye(length(xc),length(xc)); % identity matrix
A   = M + hd*alpha*B'*B;
fc  = A*xc;
u0  = zeros(size(xc));

H   = @(x) x + d2Smf.d2S(x,H.omega,H.m);  

[uCG,flag,relres,iter] = pcg(H,fc,1e-6,1e5);

testMG = norm( A\fc - uCG);
fprintf('|A\\fc-uCG|=%s\n',num2str(testMG));

%==============================================================================

%==============================================================================
FAIRmessage(sprintf('<%s> done',mfilename)); 
%==============================================================================
