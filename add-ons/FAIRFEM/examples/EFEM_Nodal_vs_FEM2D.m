% ==================================================================================
% (c) Lars Ruthotto 2012/07/26, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/lars-ruthotto.html
%
% Compares implementations of hyperelastic regularization for nodal and FEM
%
% Purpose: Reveal that alphaVolume must be chosen differently! Explore spy
%
% ==================================================================================
clear; clc;

% get nodes and triangulation
omega = [0,10,0,8]; m =[12,10]; 
xn    = reshape(getNodalGrid(omega,m),[],2);
xcc   = reshape(getCellCenteredGrid(omega,m),[],2);

% get some spline transformation and evaluate nodal and cell-centered
p = [5,6]; w = zeros([p,2]);  w(3,3,1) = 0.06; w(3,4,2) = -0.05;
yn = splineTransformation2D(w(:),xn(:),'omega',omega,'m',m+1,'p',p,'Q',[]);
yn = reshape(yn,[],2);
ycc = splineTransformation2D(w(:),xcc(:),'omega',omega,'m',m,'p',p,'Q',[]);
ycc = reshape(ycc,[],2);

% compute hyperelastic regularizer (discretized on nodal grid)
alphaVolume = 1;
reqOptn = {'alpha',1, 'alphaLength',1,'matrixFree',0};

% ATTENTION: alphaVolume must be divided by 4. A bug actually....
[Sv,dSv,d2Sv] = hyperElastic(yn(:)-xn(:),omega,m,'alphaVolume',alphaVolume/4,reqOptn{:});

Mesh = TriMesh3(omega,m);

figure(42); clf;
subplot(2,4,1);
yc = [yn; ycc];
triplot(Mesh.tri,yc(:,1),yc(:,2)); hold on; plot(yn(:,1),yn(:,2),'or');
title(sprintf('nodal grid, Sc=%1.2f',Sv));
subplot(2,4,5);
spy(d2Sv); title('non zero pattern (nodal discretization)');

Meshes = {@TriMesh1, @TriMesh2, @TriMesh3};
for type=1:3,
    Mesh = Meshes{type}(omega,m);
    
    if type==3, %  additional nodes in cell-centers
        xc = [xn;xcc]; yc = [yn; ycc];
    else % only nodal grid
        xc = xn; yc = yn;
    end
    
    [St,dSt,d2St] = hyperElasticFEM(yc(:)-xc(:),xc(:),Mesh,...
        'alphaVolume',alphaVolume,reqOptn{:});
    
    figure(42);
    subplot(2,4,1+type);
    triplot(Mesh.tri,yc(:,1),yc(:,2)); hold on; plot(yc(:,1),yc(:,2),'or');
    title(sprintf('mesh, mode=%d, Sc=%1.2f',type,St));
    subplot(2,4,5+type);
    spy(d2St); title(sprintf('non zero pattern, FEM, mode=%d',type));
end
