%==============================================================================
% This code is part of the Matlab-based toolbox 
% LagLDDDM - A Lagrangian Gauss--Newton--Krylov Solver for Mass- and 
%                        Intensity-Preserving Diffeomorphic Image Registration
% 
% For details and license info see 
% - https://github.com/C4IR/FAIR.m/tree/master/add-ons/LagLDDMM
%
% function [ z ] = solveSpectral( omega, m, tspan, nt, s1, s2, v, varargin )
%
% Solves a linear systems involving structured matrices, such as in
%
%   diffusion ST :
%
%       L = Dx'*Dx+Dy'*Dy  or  L = Dx'Dx+Dy'*Dy+Dz'*Dz ,
%
%       where Dx, Dy, Dz are partial derivative operators
%
%
%   curvatureST:
%       L = D2x*D2x+D2y*D2y or  D2x*D2x+D2y*D2y+D2z*D2z
%
%       where D2x, D2y, D2z are second-order partial derivatives
%
% For each dimension a flag determines if the corresponding term Di'*Di is
% included in the sum of L or not. A flag [1 0 1] for example will result
% in L = Dx'*Dx+Dz'*Dz.
%
% This function uses a diagonalization of L via a discrete cosine transform
% (dct) and persistent variables, only rebuilding opeRators when
% the discretization is changed (see eigLaplacian).
%
% The linear system solved has the form
% (s1*I+s2*L)z = v
%
% Input
%
%   omega  - comptuational domain
%   m      - discretization size [m1 m2] or [m1 m2 m3]
%   tspan  - time interval
%   nt     - number of cells in time
%   s1     - scaling parameter for the shift
%   s2     - scaling parameter for the Laplacian
%   v      - right hand side vector
%   varargin (see below)
%
% Output
%   z   - solution of the linear system
%==============================================================================

function [ z ] = solveSpectral( omega, m, tspan, nt, s1, s2, v, varargin )

if nargin==0
    help(mfilename)
    return
end

persistent mOld omegaOld tspanOld ntOld eigInv alphaOld regOld s1Old;

dim = numel(omega)/2;
alpha = regularizer('get','alpha');
reg  = regularizer;
dimFlag = ones(dim+1,1);
%default parameters

% overwrites default parameter
for k=1:2:length(varargin)
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;
if isempty(regularizer)
    error('%s - regularizer must be specified by FAIR or using varargin',mfilename)
end
h  = (omega(2:2:end)-omega(1:2:end))./m;
hd = prod(h);
dim = numel(omega)/2;
if not(isempty(tspan)) && not(isempty(nt))
    dt = abs(tspan(2)-tspan(1))/nt;
else
    nt = 0;
end

% restrict flags to relevant dimensions and compute effective grid size
% dimFlag = dimFlag(1:dim);
rebuild = isempty(mOld)||isempty(omegaOld)||isempty(s1Old)||...
    isempty(eigInv)||isempty(alphaOld)||isempty(regOld)...
    ||~isequal(s1,s1Old)||...
    ~isequal(m,mOld)||~isequal(omega,omegaOld)||numel(omega)~=numel(omegaOld)||...
    ~isequal(alpha,alphaOld)||~strcmp(reg,regOld) ||...
    (~isempty(tspanOld)&&~isempty(tspan)&&~isequal(tspan,tspanOld)) || ...
    (~isempty(ntOld)&&~isempty(nt)&&~isequal(nt,ntOld));
if rebuild
    fprintf('%s - rebuild inverse of structured matrix\n',mfilename);
    D = eigLaplacian(omega,m);
    switch reg
        case {'mbDiffusionCC','mfDiffusionCC','mbDiffusionST','mfDiffusionST'}
            D = alpha(1)*hd*D;
        case {'mbCurvature','mfCurvature','mbCurvatureST','mfCurvatureST'}
            D = alpha(1)*hd*D.^2;
    end

    switch reg
        case {'mbDiffusionCC','mfDiffusionCC','mbCurvature','mfCurvature'}
            eigInv = reshape(1./(s1+s2*D),m);
            eigInv(isinf(eigInv) | isnan(eigInv)) = 1;

        case {'mbDiffusionST','mfDiffusionST', 'mbCurvatureST','mfCurvatureST'}
            D = dt*D;

            eig1 = reshape(1./(s1+.5*s2*D),[],1);
            eig1(isinf(eig1) | isnan(eig1)) = 1;
            eig2 = reshape(1./(s1+s2*D),[],1);
            eig2(isinf(eig2) | isnan(eig2)) = 1;
            eigInv = [eig1(:) eig2(:)];
    end

    % save description of inverse
    mOld = m; omegaOld = omega; tspanOld=tspan; ntOld = nt; alphaOld=alpha;
    regOld=reg; s1Old = s1;
end


switch reg
    case {'mbDiffusionCC','mfDiffusionCC','mbCurvature','mfCurvature'}
        getEig = @(k) eigInv;
    case {'mbDiffusionST','mfDiffusionST', 'mbCurvatureST','mfCurvatureST'}
        getEig = @(k) reshape(eigInv(:,1+((k>1)&&(k<nt+1))),m);
end
% discrete cosine transform diagonalized inversion
nrhs = size(v,2);
if dim == 3
    if nrhs==1
        V = reshape(full(v),[m dim nt+1]);
        Z = zeros([m dim nt+1]);
        for d=1:dim
            for k=1:nt+1
                Z(:,:,:,d,k) = idctn(getEig(k).*dctn(squeeze(V(:,:,:,d,k)),'dimFlag',dimFlag),'dimFlag',dimFlag);
            end
        end
        z = Z(:);
    else
        error('nyi')
    end
else
    if nrhs==1
        V = reshape(full(v), [m dim nt+1]);
        Z = zeros([m dim nt+1]);
        for d=1:dim
            for k=1:nt+1
                Z(:,:,d,k) = idctn(getEig(k).*dctn(squeeze(V(:,:,d,k)),'dimFlag',dimFlag),'dimFlag',dimFlag);
            end
        end
        z = Z(:);
    else
        error('nyi')
    end
end



