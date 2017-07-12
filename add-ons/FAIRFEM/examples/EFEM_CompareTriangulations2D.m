% ==================================================================================
% (c) Lars Ruthotto 2012/07/26, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/lars-ruthotto.html
%
% 2D academic example for hyperelastic image registration using a finite element method (FEM)
%
% Purpose: Compare how different triangulations vary the registration result
%
% Data: This academic example is based on Gary Christensens famous Disc to C example.
%
% Q: Is data symmetric
%    are coefficients symmetric
%    Wo bricht symmetrie?
%
% ==================================================================================

setup2Ddisc2CData

% initialize modules
imgModel('reset','imgModel','splineInter','regularizer','moments','theta',1e-2);
aVol = 1;
regularizer('reset','regularizer','mbHyperElasticFEM','alpha',2e2,'alphaVolume',aVol);
level =4; omega = ML{level}.omega; m = ML{level}.m;
[T,R] = imgModel('coefficients',ML{level}.T,ML{level}.R,omega(1,:),'out',0);

n2s = @(n) num2str(n);
% enforce symmetry of data
R(9:end,:) = flipud(R(1:8,:));
T(9:end,:) = flipud(T(1:8,:));

%% run FEM implementation
Meshes = {@TriMesh1,@TriMesh2,@TriMesh3};
for type=1:3, % run over different triangulations
    Mesh = Meshes{type}(omega,m);
    xn   = Mesh.xn; tri = Mesh.tri;
    
    % interpolate reference image in centers of triangles
    Rc = imgModel(R,omega,Mesh.mfPi(xn,'C'));
    
    % open FAIRplotsFEM
    FAIRplotsFEM('reset','mode','FEM','fig',level,'plots',1);
    FAIRplotsFEM('init',struct('Tc',T,'Rc',R,'Mesh',Mesh));
    
    % setup objective function
    NPIRfctn = @(yc) FEMobjFctn(T,Rc,Mesh,xn,yc);
    % report status
    NPIRfctn([]);
    % optimize
    [yc,his] = GaussNewton(NPIRfctn,xn(:),'yc',xn(:),'yStop',xn(:),'maxIter',80,...
        'Plots',@FAIRplotsFEM,'lineSearch',@ArmijoDiffeomorphicFEM);
    
    % store results
    eval(['MeshType' n2s(type) '=Mesh;yType' n2s(type) '=yc;']);
end

%% run NPIR version
regularizer('set','regularizer','mbHyperElastic','alphaVolume',aVol/4);
xc = getNodalGrid(omega,m);
Rc = imgModel(R,omega,center(xc,m));
% setup FAIR plots
FAIRplots('reset','mode','NPIR','fig',level);
FAIRplots('init',struct('Tc',T,'Rc',R,'omega',omega,'m',m));
% setup objective function
NPIRfctn = @(yc) NPIRobjFctn(T,Rc,omega,m,xc,yc);
% report status
NPIRfctn([]);
% optimize
[yNPIR,his] = GaussNewton(NPIRfctn,xc,'yc',xc,'yStop',xc,'maxIter',80,...
    'Plots',@FAIRplots,'lineSearch',@ArmijoDiffeomorphic);
% prepare results
ytriNPIR = [reshape(yNPIR,[],2);reshape(nodal2center(yNPIR,m),[],2)];
%% display results
clc;
fprintf('Evaluating Mesh Quality\n\n%6s %10s %10s %10s %10s %10s %10s\n',...
    'Mesh','min(q1)','mean(q1)', 'min(q2)','mean(q2)','min(q3)','mean(q3)');
figure(42); clf;
for type=1:3,
    eval(['Mesh=MeshType' n2s(type) ';yc=yType' n2s(type) ';']);
    tri = Mesh.tri; xn = Mesh.xn;
    % display triangulation
    subplot(2,4,type);
    triplot(tri,xn(:,1),xn(:,2)); hold on; plot(xn(:,1),xn(:,2),'or');
    axis(omega); title(sprintf('xn, mode=%d',type));
    
    subplot(2,4,4+type);
    yc = reshape(yc,[],2);
    viewImage2Dsc(dataT,omega,size(dataT));
    colormap gray
    C = colormap; colormap(flipud(C));
    hold on;
    triplot(tri,yc(:,1),yc(:,2),'linewidth',3); hold on; plot(yc(:,1),yc(:,2),'.r','MarkerSize',20);
    axis(omega); title(sprintf('yc, mode=%d',type));
%     [q1 q2 q3] = getMeshQuality(Mesh,yc);
%     fprintf('%6s %10s %10s %10s %10s %10s %10s\n',...
%         n2s(type),n2s(min(q1)),n2s(mean(q1)), n2s(min(q2)),n2s(mean(q2)),...
%         n2s(min(q3)),n2s(mean(q3)));
%     
end
% add NPIR results
subplot(2,4,4);
triplot(MeshType3.tri,MeshType3.xn(:,1),MeshType3.xn(:,2)); hold on; 
plot(MeshType1.xn(:,1),MeshType1.xn(:,2),'or');
axis(omega); title(sprintf('xn, NPIR'));

subplot(2,4,8);
triplot(MeshType3.tri,ytriNPIR(:,1),ytriNPIR(:,2));
yNPIR = reshape(yNPIR,[],2);
hold on; plot(yNPIR(:,1),yNPIR(:,2),'or');
axis(omega); title(sprintf('yc, NPIR'));
% [q1 q2 q3] = getMeshQuality(MeshType3,ytriNPIR);
% fprintf('%6s %10s %10s %10s %10s %10s %10s\n',...
%     'NPIR',n2s(min(q1)),n2s(mean(q1)), n2s(min(q2)),n2s(mean(q2)),...
%     n2s(min(q3)),n2s(mean(q3)));


return;
%% print plots for thesis
close all;
%% show initial data
Tc = imgModel(T,omega,getCellCenteredGrid(omega,m));
Rc = imgModel(R,omega,getCellCenteredGrid(omega,m));

viewImg = @(I) viewImage2Dsc(I,omega,m,'colormap',flipud(gray(256)));
figure(1); clf;
viewImg(Tc);
figure(2); clf;
viewImg(Rc);

%% show transformations and interpolated templates
for type=1:3,
    eval(['Mesh=MeshType' n2s(type) ';yc=yType' n2s(type) ';']);
    tri = Mesh.tri; xn = Mesh.xn;
    yc = reshape(yc,[],2);
    
    figure(2+type);clf;
    viewImg(Tc); hold on;
    
    tp= triplot(tri,yc(:,1),yc(:,2),'LineWidth',2,'Color','y'); hold on; 
    plot(yc(:,1),yc(:,2),'ob','LineWidth',2);
    axis(omega);   
end


%% show interpolated templates
for type=1:3,
    eval(['Mesh=MeshType' n2s(type) ';yc=yType' n2s(type) ';']);
    tri = Mesh.tri; xn = Mesh.xn;
    yc = reshape(yc,[],2);
    
    figure(5+type);clf;
    viewImg(Tc); hold on;
    Topt = imgModel(T,omega,Mesh.mfPC(yc));
    ph = trisurf(tri,xn(:,1),xn(:,2),0*yc(:,2),Topt,'EdgeColor','k','LineWidth',.5);
    view([0 90])    ; axis(omega);axis xy;
end

%% print images
outdir = '/Users/larsruthotto/Dropbox/2012-Diss/plots/EFEM_CompareTriangulations';
names = {'dataT','dataR','yType1','yType2','yType3','ToptType1','ToptType2','ToptType3'};
for i=3:size(names,2),
   figure(i);
   set(figure(i),'PaperPositionMode','auto');
   set(gca,'Position',[0 0 1 1]);
   set(gcf,'Units','inches');
   set(gcf,'Position',[0 0 10 10]);
   axis off;
   printFigure(figure(i),fullfile(outdir,['CompareTriangulations-' names{i} '.pdf']),...
       'printOpts','-dpdf','printFormat','.pdf');
    
end

