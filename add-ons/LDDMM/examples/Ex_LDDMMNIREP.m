close all; clc; clear all;

setupNIREPDataNA02;
%setupNIREPDataNA03;
%setupNIREPDataNA04;
%setupNIREPDataNA05;
%setupNIREPDataNA06;

imgModel('reset','imgModel','linearInterMex');
trafo('reset','trafo','rigid3D');
distance('reset','distance','SSD');
viewImage('reset','viewImage','viewOrthoSlices2D');

profilecode = false;
% example specific parameters
switch example
    case ['3D-nirep-',dataset]
        lvl = 4;
        pad = 1;
        mV = @(m) m;
        nt = 1;
        ntpost = 20;
        N  = 5;  % default was 3
        minLevel = 5;
        maxLevel = 7;
        maxIterNPIR = 50;
    otherwise
        error('MP-LDDMM example %s - nyi');
end

tolJ=5E-3;
%shiftlist  = [0, 1E-8, 1E-6, 1E-4, 1E-2];
shiftlist  = [1E-2];
alpha1 = [5,10,25,50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000];
alpha2 = 0;

for j = 1:numel(shiftlist)
    shift = shiftlist(j);
    for i = 5:9
        alpha = [alpha1(i),alpha2];
        % 1) setup grid for velocities (padded)
        omegaV = omega; omegaV(1:2:end) = omegaV(1:2:end)-pad;  omegaV(2:2:end) = omega(2:2:end)+pad;

        % 2) setup regularizer (and decide for stationary or nonstationary velocity)
        regularizer('reset','regularizer','mfDiffusionST','alpha',alpha,'nt',nt,'HessianShift',shift,'fig',0); % nonstationary velocity
%        regularizer('reset','regularizer','mfDiffusionCC','alpha',alpha(1),'HessianShift',shift,'fig',0); % stationary velocity

%        regularizer('reset','regularizer','mfCurvatureST','alpha',alpha,'nt',nt,'HessianShift',shift,'fig',0) % nonstationary velocity
%        regularizer('reset','regularizer','mfCurvature','alpha',alpha(1),'HessianShift',shift,'fig',0); % stationary velocity

        % check if we solve stationary
        if (   strcmp(regularizer,'mfDiffusionCC') ...
             | strcmp(regularizer,'mfCurvature') )
            nt = 0;
        end

        solve = true;
        postproc = true;

        %% construct filename
        if strcmp(regularizer,'mfHyperElastic')
            filename = ['HyperElastic-',example,'-',regularizer,'-alpha-',num2str(alpha(1)),'-',num2str(alpha(2)),'-',num2str(alpha(3)),'-',num2str(alpha(4))];
        else
            filename = ['LDDMM-',example,'-',regularizer,'-alpha-',num2str(alpha(1)),'-',num2str(alpha(2)),'-N-',num2str(N),'-nt-',num2str(nt),'-shift-',num2str(shift),'-tolJ-',num2str(tolJ)];
        end

        if exist(['results/',filename,'.mat'], 'file')
            solve = false;
            postproc = false;
        end


        if (solve)
            disp(['solving ',filename]);

            if (profilecode)
                profile on;
            end
            diary(fullfile('results',[filename,'.log']));
            [vc,yc,wc,his,para] = MLLDDMM(MLdata,'minLevel',minLevel,'maxLevel',maxLevel,'omegaV',omegaV,...
                'mV',mV,'N',N,'parametric',0,'maxIterNPIR',maxIterNPIR,'tolJ',tolJ,'solverNPIR',[],'fig',0);
            diary('off')

            if (profilecode)
                profile off;
                profsave(profile('info'),'lddmmprofile_results')
            end
            [~,interOpts] = inter;
            [~,distOpts ] = distance;
            [~,regOpts  ] = regularizer;

            ntt = round(numel(vc)/(3*prod(m)))-1;
            if (nt ~= ntt), error('nt is wrong'); end
            xc = getCellCenteredGrid(omega,m);
            if nt < 1
                yc = getTrafoFromVelocityRK4(vc,xc,'omega',omegaV,'m',mV(MLdata{maxLevel}.m),'N',N,'tspan',[1 0]);
                yInv = getTrafoFromVelocityRK4(vc,xc,'omega',omegaV,'m',mV(MLdata{maxLevel}.m),'N',N,'tspan',[0 1]);
            else
                yc = getTrafoFromInstationaryVelocityRK4(vc,xc,'omega',omegaV,'m',mV(MLdata{maxLevel}.m),'N',N,'nt',nt,'tspan',[1 0]);
                yInv = getTrafoFromInstationaryVelocityRK4(vc,xc,'omega',omegaV,'m',mV(MLdata{maxLevel}.m),'N',N,'nt',nt,'tspan',[0 1]);
            end
            save(fullfile('results',[filename,'.mat']),'vc','yc',...
                'wc','his','mV','minLevel','maxLevel','para','interOpts','distOpts','regOpts','N','yInv')
        end

        if postproc
            % generate figures
            close all;
            load(fullfile('results',[filename,'.mat']));
            dim  = numel(omega)/2;

            xc = getCellCenteredGrid(omega,m);
            if nt < 1
                yc = getTrafoFromVelocityRK4(vc,xc,'omega',omegaV,'m',mV(MLdata{maxLevel}.m),'N',ntpost,'tspan',[1 0]);
                yInv = getTrafoFromVelocityRK4(vc,xc,'omega',omegaV,'m',mV(MLdata{maxLevel}.m),'N',ntpost,'tspan',[0 1]);
            else
                yc = getTrafoFromInstationaryVelocityRK4(vc,xc,'omega',omegaV,'m',mV(MLdata{maxLevel}.m),'N',ntpost,'nt',nt,'tspan',[1 0]);
                yInv = getTrafoFromInstationaryVelocityRK4(vc,xc,'omega',omegaV,'m',mV(MLdata{maxLevel}.m),'N',ntpost,'nt',nt,'tspan',[0 1]);
            end

            computeOverlap;

            figDir = ['results/',filename] ;
            if ~exist(figDir,'dir'), mkdir(figDir); end;
            dim  = numel(omega)/2;

            numfig = 7;
            figid = 1;
            viewImage('reset','viewImage','viewOrthoSlices2D','colormap','gray');
            fig = figure(figid); clf;
            fig.Name = 'dataR';
            viewImage(dataR,omega,m);
            cax = caxis;

            figid = figid + 1;
            fig = figure(figid); clf;
            fig.Name = 'dataT';
            viewImage(dataT,omega,m);
            caxis(cax);

            figid = figid + 1;
            fig = figure(figid); clf;
            fig.Name = 'res0';
            viewOrthoSlices2D(abs(dataR-dataT),omega,m,'color1',[0 0 0],'color2',[0 0 0]);
            colormap gray;
            colormap(flipud(colormap));
            caxd = cax;

            Topt = imgModel(dataT,omega,center(yc,m));
            figid = figid + 1;
            fig = figure(figid); clf;
            fig.Name = 'Topt';
            viewImage(Topt,omega,m);
            caxis(cax);

            figid = figid + 1;
            fig = figure(figid); clf;
            fig.Name = 'resOpt';
            viewOrthoSlices2D(abs(dataR(:)-Topt),omega,m,'color1',[0 0 0],'color2',[0 0 0]);
            caxis(caxd);
            colormap gray;
            colormap(flipud(colormap));

            figid = figid + 1;
            fig = figure(figid); clf;
            fig.Name = 'vol';
            Jac = geometry(yc,m-1,'Jac','matrixFree',1,'omega',omega);
            viewOrthoSlices2D(Jac,omega,m-1,'color1',[0 0 0],'color2',[0 0 0]);
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

            figid = figid + 1;
            fig = figure(figid); clf;
            fig.Name = 'volInv';
            JacInv = geometry(yInv,m-1,'Jac','matrixFree',1,'omega',omega);
            viewOrthoSlices2D(JacInv,omega,m-1,'color1',[0 0 0],'color2',[0 0 0]);
            colormap jet;
            cb = colorbar('East');
            cpos = cb.Position;
            cb.Position(4) = cpos(4)*.4;
            cb.Position(2) = .56;
            cb.Ticks = [min(JacInv) max(JacInv)];
            caxis([min(JacInv) max(JacInv)]);
            cb.TickLabels = {sprintf('%1.2f',min(JacInv)), sprintf('%1.2f',max(JacInv))};
            cb.FontSize = 50;
            cb.Color = 'k';

            if (figid ~= numfig), error('figure setup is wrong'); end;

            str = {'dataR','dataT','res0','Topt','resOpt','vol','volInv','volErr'};
            for k=1:numfig
                figure(k);
                axis off;
                title([]);
                printFigure(gcf,fullfile(figDir,['LDDMM-' str{k} '_' regularizer '.png']),...
                    'printOpts','-dbmp','printFormat','.bmp');
            end

            for k=1:nt+1
                fig = figure(9+k); clf;
                fig.Name = sprintf('vel(%d)',k);
                vcc = reshape(vc,[m,3,nt+1]);
                vcc = reshape(squeeze(vcc(:,:,end/2,1:2,:)),[],nt+1);
                viewVelocity(vcc(:,k),[omega(1),omega(2),omega(3),omega(4)],[m(1),m(2)]);
            end
            for k=1:nt+1
                figure(9+k);
                axis off;
                title([]);
                printFigure(gcf,fullfile(figDir,['LDDMM_vc-' num2str(k) '_' regularizer '.png']),...
                    'printOpts','-dbmp','printFormat','.bmp');
            end
        end
    end
end
