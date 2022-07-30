%%Traces from example subject to show how data is summarized
%% Load data
% load('/Volumes/Users/Dulce/R01_Nimbus2021/VROG_Devon/VrG_Devon.mat')
% subID = 'CTR_01';
% scriptDir = fileparts(matlab.desktop.editor.getActiveFilename); 
% load([scriptDir '/data/' subID])

%% Set muscle to plot



%% Calf muscle 

normalize = 1;  % 1 to normalize data
normCond = {'TM mid 1'};
conds={'TM mid','Adaptation','Adaptation','Pos Short','Neg Short'};
late=[1  0 1 1 0 0];
strides=[40 40 40 20 20 ];
IgnoreStridesEarly=[1  50 1 50 50];
% fh=figure(1)
muscle={'TA', 'MG'};
plotEMGtraces_CI(expData,conds,muscle,late,strides,normalize,normCond,IgnoreStridesEarly);

muscle={'TA', 'LG'};
plotEMGtraces_CI(expData,conds,muscle,late,strides,normalize,normCond,IgnoreStridesEarly);

muscle={'TA', 'PER'};
plotEMGtraces_CI(expData,conds,muscle,late,strides,normalize,normCond,IgnoreStridesEarly);

%% Thig muscles 
muscle={'VL', 'BF'};
normalize = 1;  % 1 to normalize data
normCond = {'TM mid 1'};
conds={'TM mid 1','Adaptation','Adaptation'};
late=[1 0 1];
strides=[40 40 40];
IgnoreStridesEarly=[1 50 1];
% fh=figure(1)
plotEMGtraces_CI(expData,conds,muscle,late,strides,normalize,normCond,IgnoreStridesEarly);
