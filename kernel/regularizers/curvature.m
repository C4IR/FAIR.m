%==============================================================================
% This code is part of the Matlab-based toolbox
% FAIR - Flexible Algorithms for Image Registration. 
% For details see 
% - https://github.com/C4IR and
% - http://www.siam.org/books/fa06/
%==============================================================================
%
% function [Sc,dS,d2S] = curvature(uc,omega,m,varargin)
%
% computes curvature regularization energy for uc = yc - yRef
% (cell centered)
%
% S(u) = 0.5 * \int_{\omega} u(x)' * B' * B * u(x) dx,
%
% where B is the curvature operator, see getCurvatureMatrix .
%
% Input:
% ------
%   uc           displacement field (staggered)
%   omega       spatial domain
%   m             number of discretization points
%   varargin    optional parameters (see below)
%
% Output:
% -------
%   Sc          current value  (0.5 * hd * uc'*B'*B*uc)
%   dS          derivative     (hd * uc'*B'*B)
%   d2S           Hessian        (B'*B)
%  if ~matrixFree,  d2S is sparse matrix; else, d2S is struct endif        
%
% =============================================================================

function [Sc,dS,d2S] = curvature(uc,omega,m,varargin)
if nargin == 0
    help(mfilename);
    runMinimalExample;
    Sc = 'endOfMinimalExample';
    return;
end

persistent A omegaOld mOld alphaOld

% if ~exist('mOld','var'),     mOld = [];     end;
% if ~exist('omegaOld','var'), omegaOld = []; end;
% if ~exist('alphaOld','var'), alphaOld = []; end;

matrixFree  = 0;
alpha       = 1;

for k=1:2:length(varargin), % overwrites default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

hd  = prod((omega(2:2:end)-omega(1:2:end))./m);

if not(matrixFree), % matrix-based
  build = isempty(mOld) || isempty(omegaOld) || isempty(alphaOld) ...
    || length(mOld) ~= length(m) || length(omegaOld) ~= length(omega) ...
    || any(mOld ~= m) || any(omegaOld~=omega) || (alpha ~=alphaOld);
  
  if build,
    mOld = m; omegaOld = omega; alphaOld = alpha;
    A = getCurvatureMatrix(omega,m);
    A = alpha*hd*(A'*A);
  end;
  dS  = uc'*A;
  Sc  = 0.5*dS*uc;
  d2S = A;
else % matrix-free
  d2S.regularizer = regularizer;
  d2S.alpha  = alpha;
  d2S.By     = @(u,omega,m) curvatureOperator(u,omega,m,'By');
  d2S.BTy    = @(u,omega,m) curvatureOperator(u,omega,m,'BTy');
  d2S.B      = @(omega,m)   getCurvatureMatrix(omega,m);
  
  if (length(m) == 2) || (length(m) == 3) %mex solution is available
    [Sc, dS] = curvatureMexC(uc,omega,m,alpha);
    d2S.diag = @(omega,m) curvatureDiagMex(omega,m,alpha);
    d2S.d2S  = @(uc,omega,m) alpha * hd * curvatureHessianMex(uc,omega,m);
  else
    d2S.d2S  = @(uc,omega,m) ...
      alpha * prod((omega(2:2:end)-omega(1:2:end))./m) * ...
      curvatureOperator(curvatureOperator(uc,omega,m,'By'),omega,m,'BTy');
    
    d2S.diag = @(omega,m) getDiag(omega,m,alpha);
    dS   = d2S.d2S(uc,omega,m)';
    Sc   = .5*dS*uc;
  end
    
end

%------------------------------------------------------------------------------

function D = getDiag(omega,m,alpha)
%Compute Diagonal of alpha*hd*B'*B = d2S
h    = (omega(2:2:end)-omega(1:2:end))./m;
hd   = prod(h);

B=getCurvatureMatrixSmall(omega,m);

D=alpha*hd*sum(B.*B,2);
D=repmat(D,dim,1);

%------------------------------------------------------------------------------

function By = curvatureOperator(uc,omega,m,flag)
dim = length(omega)/2;
flag = sprintf('%s-%dD',flag,dim);
switch flag,
  case {'By-2D','BTy-2D'},
    %
    % By = | \Delta y^1 0          | = | d_1^2 y^1 + d_2^2 y^1 |
    %      | 0          \Delta y^2 | = | d_1^2 y^2 + d_2^2 y^2 |   
    uc = reshape(uc,[m,2]);
    By = zeros(size(uc));
    d2 = @(i) D2(i,omega,m);
    By(:,:,1) = d2(1)*uc(:,:,1) + uc(:,:,1)*d2(2);
    By(:,:,2) = d2(1)*uc(:,:,2) + uc(:,:,2)*d2(2);
    By = By(:);
  
  case {'By-3D','BTy-3D'},
    %
    %      | \Delta y^1 0          0          | = | d_1^2 y^1 + d_2^2 y^1 + d_3^2 y^1 |
    % By = | 0          \Delta y^2 0          | = | d_1^2 y^2 + d_2^2 y^2 + d_3^2 y^2 |   
    %      | 0          0          \Delta y^3 | = | d_1^2 y^3 + d_2^2 y^3 + d_3^2 y^3 |   

    % efficient implementation of Bcurvature*Y    
    n  = prod(m);            % number of voxels
    uc = reshape(uc,[m,3]);  % note that now uc(:,:,:,ell) is the ell-th
                             % component of Y=(Y^1,Y^2,Y^3)
                             % of size m(1)-by-m(2)-by-m(3)
    By = zeros(numel(uc),1); % allocate memory for the output
    d2 = @(i) D2(i,omega,m); % this is a shortcut to the discrete second derivative

    % the following line is a shortcut for 
    %  - permuting the 3d-array using the permutation J
    %  - reshape it to a 2D-array of size q-by-prod(m)/q, where q=m(J(1))
    %  - multiply by A (which is q-by-q)
    %  - undo the reshape, i.e. make the result to m(J(1))-by-m(J(2))-by-m(J(3))
    %  - undo the permutation
    operate = @(A,z,J) ipermute(reshape(A*reshape(permute(z,J),m(J(1)),[]),m(J)),J);

    % run over all compunents y^ell of Y=(Y^1,Y^2,Y^3)
    for ell=1:3,
      % compute
      % (I_3\otimes I_2\otimes d2(1) + I_3\otimes d2(2)\otimes I_1 ...
      % + d2(3)\otimes I_2\otimes I_1) y^ell
      for k=1:3,
        z = operate(d2(k),uc(:,:,:,ell),[k,setdiff(1:3,k)]);
        By((ell-1)*n+(1:n)) = By((ell-1)*n+(1:n)) + reshape(z,[],1);
      end;
    end;   
end;

%------------------------------------------------------------------------------

function D = D2(i,omega,m)
h = (omega(2:2:end)-omega(1:2:end))./m;
D = spdiags(ones(m(i),1)*[1,-2,1],-1:1,m(i),m(i))/h(i)^2;
D([1,end]) = -D([2,end-1]);

%------------------------------------------------------------------------------

function runMinimalExample
  % --- 2D --- 
  omega = [0 2 0 2];
  m     = [17 19];
  xc    = getCellCenteredGrid(omega,m);
  uc    = 1e-2*randn(size(xc));
  [mbS, mbdS, mbd2S] = feval(mfilename,uc,omega,m,'alpha',1,'matrixFree',0);
  [mfS, mfdS, mfd2S] = feval(mfilename,uc,omega,m,'alpha',1,'matrixFree',1);
  fprintf('Curvature regularizer : %f (mb) | %f (mf) \n' , mbS,mfS);
  RE = norm(mbdS-mfdS)/norm(mbdS);
  fprintf('rel.error dS (mb vs. mf) : %e \n' , RE);
  figure(1)
  subplot(1,2,1);
  spy(mbd2S)
  title(sprintf('%s-spy(d2S)-(2D)',mfilename));
  
  % --- 3D --- 
  omega = [0 2 0 3 0 2];
  m     = [6 5 7];
  xc    = getCellCenteredGrid(omega,m);
  uc    = 1e-2*randn(size(xc));
  [mbS, mbdS, mbd2S] = feval(mfilename,uc,omega,m,'alpha',1,'matrixFree',0);
  [mfS, mfdS, mfd2S] = feval(mfilename,uc,omega,m,'alpha',1,'matrixFree',1);
  fprintf('Curvature regularizerl : %f (mb) | %f (mf) \n' , mbS,mfS);
  RE = norm(mbdS-mfdS)/norm(mbdS);
  fprintf('rel.error dS (mb vs. mf) : %e \n' , RE);
  figure(1)
  subplot(1,2,2);
  spy(mbd2S)
  title(sprintf('%s-spy(d2S)-(3D)',mfilename));
% =============================================================================


  