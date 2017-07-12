close all; clc; clear all;

setupNIREPDataNA02;
%setupNIREPDataNA03;
%setupNIREPDataNA04;
%setupNIREPDataNA05;

imgModel('reset','imgModel','linearInterMex');
trafo('reset','trafo','rigid3D');
distance('reset','distance','SSD');
viewImage('reset','viewImage','viewOrthoSlices2D');

% example specific parameters
switch example
    case ['3D-nirep-',dataset]
        alpha(1) = 100;
        alpha(2) = 1;
        alpha(3) = .1;
        alpha(4) = 1;
        pad = 0;
        minLevel = 5;
        maxLevel = 7;
        maxIterNPIR = 50;
    otherwise
        error('MP-LDDMM example %s - nyi');
end

Dshow = @(T,R,omega,m) viewIP(abs(T-R),omega,m,'colormap',gray(256));
FAIRplots('set','Dshow',Dshow);


solve = true;
postproc = true;



% 2) setup regularizer
regularizer('reset','regularizer','mfHyperElastic','alpha',alpha(1),...
    'alphaLength',alpha(2),'alphaArea',alpha(3),'alphaVolume',alpha(4));

filename = ['HyperElastic-',example,'-',regularizer,'-alpha-',num2str(alpha(1)),'-',num2str(alpha(2)),'-',num2str(alpha(3)),'-',num2str(alpha(4))];

if (solve)
    diary(fullfile('results',[filename,'.log']));
    % finally: run the MultiLevel Non-Parametric Image Registration
    [yc,wc,his] = MLIR(MLdata,'parametric',0,'NPIRLS',@ArmijoBacktrack,'minLevel',minLevel,'maxLevel',maxLevel,'maxIterNPIR',maxIterNPIR);
    diary('off')
    [~,interOpts] = inter;
    [~,distOpts ] = distance;
    [~,regOpts  ] = regularizer;

    save(fullfile('results',[filename,'.mat']),'yc','wc','his','m','minLevel','maxLevel','interOpts','distOpts','regOpts')
end

if postproc
    % generate figures
    close all;
    load(fullfile('results',[filename,'.mat']));

    computeOverlap;

    % generate figures
    figDir = ['results/',filename];
    if ~exist(figDir,'dir'), mkdir(figDir); end;
    dim  = numel(omega)/2;

    viewImage('reset','viewImage','viewOrthoSlices2D','colormap','gray');
    fig = figure(1); clf;
    fig.Name = 'dataR';
    viewImage(dataR,omega,m);
    cax = caxis;

    fig = figure(2); clf;
    fig.Name = 'dataT';
    viewImage(dataT,omega,m);
    caxis(cax);

    fig = figure(3); clf;
    fig.Name = 'res0';
    % viewImage(abs(dataR-dataT),omega,m)
    viewOrthoSlices2D(abs(dataR-dataT),omega,m,'color1',[0 0 0],'color2',[0 0 0]);
    colormap gray;
    colormap(flipud(colormap));
    % caxd = caxis;
    caxd = cax;

    Topt = imgModel(dataT,omega,center(yc,m));
    fig = figure(4); clf;
    fig.Name = 'Topt';
    viewImage(Topt,omega,m);
    % colormap gray
    caxis(cax);

    fig = figure(5); clf;
    fig.Name = 'resOpt';
    % viewImage(abs(dataR(:)-Topt),omega,m)
    viewOrthoSlices2D(abs(dataR(:)-Topt),omega,m,'color1',[0 0 0],'color2',[0 0 0]);
    caxis(caxd);
    colormap gray;
    colormap(flipud(colormap));

    fig = figure(6); clf;
    fig.Name = 'vol';
    Jac = geometry(yc,m,'Jac','matrixFree',1,'omega',omega);
    viewOrthoSlices2D(Jac,omega,m,'color1',[0 0 0],'color2',[0 0 0]);
    colormap jet;
    cb = colorbar('East');
    cpos = cb.Position;
    cb.Position(4) = cpos(4)*.4;
    cb.Position(2) = .56;
    cb.Ticks = [min(Jac) max(Jac)];
    caxis([min(Jac) max(Jac)]);
    cb.TickLabels = {sprintf('%1.2f',min(Jac)), sprintf('%1.2f',max(Jac))};
    cb.FontSize = 50;
    cb.Color = 'k';

    %str = {'dataR','dataT','res0','Topt','resOpt','yOpt','yInv','vol','T+yc'};
    str = {'dataR','dataT','res0','Topt','resOpt','vol'};
    for k=1:6
        figure(k);
        axis off;
        title([]);
        printFigure(gcf,fullfile(figDir,['LDDMM-' str{k} '_' regularizer '.png']),...
            'printOpts','-dbmp','printFormat','.bmp');
    end
end
