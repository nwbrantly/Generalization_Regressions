%%Traces from example subject to show how data is summarized
%% Load data
% load('/Volumes/Users/Dulce/R01_Nimbus2021/VROG_Devon/VrG_Devon.mat')
% subID = 'CTR_01';
% scriptDir = fileparts(matlab.desktop.editor.getActiveFilename); 
% load([scriptDir '/data/' subID])

%% Set muscle to plot

muscle={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT', 'VM', 'VL', 'RF', 'TFL', 'GLU','HIP'};
normalize = 1;  % 1 to normalize data
normCond = {'TM base'};

%% Baseline condtions 
conds={'OG base','TM fast','TM base','TM slow'};
late=1;
strides=40;
fh=plotEMGtraces(expData,conds,muscle,late,strides,normalize,normCond);

%% Late condition 
late=1;
strides=40;
conds={'OG base','TM fast','TM slow'...
    'TM base','Adaptation',...
    'Post 1','Post 2'};
fh=plotEMGtraces(expData,conds,muscle,late,strides,normalize,normCond);


%% Early conditions 
late=0;
strides=30;
conds={'TM base','Pos short',...
    'Neg Short','Adaptation',...
    'Post 1','Post 2'};
fh=plotEMGtraces(expData,conds,muscle,late,strides,normalize,normCond);

%% Early conditions 
late=0;
strides=30;
conds={'Pos short',...
    'Neg Short','Adaptation',...
        };
fh=plotEMGtraces(expData,conds,muscle,late,strides,normalize,normCond);

%% save figures
% if late
%     if baseOnly
%         saveas(fh, [scriptDir '/EMGTraces/' subID '_BaseLate.png']);
%     else
%         saveas(fh, [scriptDir '/EMGTraces/' subID '_Late.png']);
%     end
% else
%     saveas(fh, [scriptDir '/EMGTraces/' subID '_Early.png']);
% end
