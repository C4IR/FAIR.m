%==============================================================================
% FAIR: Flexible Algorithms for Image Registration toolbox
% Jan Modersitzki, Lars Ruthotto and Fabian Gigengack
% see https://github.com/C4IR/FAIR.m for details and license
%==============================================================================
% 
% setup toolbox path
%==============================================================================

function FAIRstartup

FAIRpath = fileparts(which(mfilename('fullpath')));

% initial message
message = @(str) fprintf('%s\n',str{:});

message({
    '% =============================================================================='
    'FAIR: Flexible Algorithms for Image Registration'
    '(c) Jan Modersitzki, Lars Ruthotto'
     sprintf('Set path on [%s], FAIRpath is [%s]\n',computer,FAIRpath);
    });

% create FAIRpath, for this machine
cdAct = pwd;
addpath(cdAct);

cd(FAIRpath);
str = sprintf('function value=FAIRpath; value=''%s'';',FAIRpath);
[f,msg] = fopen(fullfile(FAIRpath,'kernel','tools','FAIRpath.m'),'w');
fprintf(f,'%s\n',str);
fclose(f);

% create temp folder for temporary data, testing and debugging, if necessary
fairTemp = fullfile(FAIRpath,'temp');
if ~exist(fairTemp,'dir'), mkdir(fairTemp);  end;

% debit = contents
% FAIRtoolboxes = debit(find(cellfun(@exist,debit)==7));

% add subfolder containing the different parts
addpath(genpath(fullfile(FAIRpath,'kernel')));
addpath(fullfile(FAIRpath,'temp'));
addpath(genpath(fullfile(FAIRpath,'add-ons')));
message({
    '[new to FAIR? try the numerous examples in kernel/examples]'
    '% =============================================================================='
    ''
    });

cd(cdAct)
%==============================================================================
