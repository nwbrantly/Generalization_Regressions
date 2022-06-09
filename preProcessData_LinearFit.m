%%
%Running this script requires labTools (github.com/pittSMLlab/labTools/)

% addpath(genpath('../../../EMG-LTI-SSM/'))
% addpath(genpath('../../../matlab-linsys/'))
% addpath(genpath('../../../robustCov/'))
%% Aux vars:
% groupName='controls';
% groupName='YoungAbuptTM';
% matDataDir='./';
% loadName=[matDataDir groupName];
% load(loadName)

%%
% group=controls;
% group=TMFullAbrupt;
% group={adaptData};
%%
clear; close all; clc;

% set script parameters, SHOULD CHANGE/CHECK THIS EVERY TIME.
groupID = 'ATR';
saveResAndFigure = false;
plotAllEpoch = true;
plotIndSubjects = true;
plotGroup = true;
bootstrap=true;

scriptDir = cd;% fileparts(matlab.desktop.editor.getActiveFilename);
files = dir ([groupID '*params.mat']);



sub = {};
subID = {};

% for i = 1:size(files,1)
% 
%         sub{end+1} = files(i).name;
%         subID{end+1} = sub{end}(1:end-10);
% 
% end
% % n_subjects = size(files,1) - session2_n_subjects;
% subID
session2_n_subjects = 0;
sub = {};
subID = {};
session2subID = {};
session2sub = {};
for i = 1:size(files,1)
    if contains(files(i).name,'Session2')
        session2_n_subjects = session2_n_subjects + 1;
        session2sub{end+1} = files(i).name;
        session2subID{end+1} = session2sub{end}(1:end-10);
    else
        sub{end+1} = files(i).name;
        subID{end+1} = sub{end}(1:end-10);
    end
end
n_subjects = size(files,1) - session2_n_subjects;
subID
session2subID

group=adaptationData.createGroupAdaptData(sub); %loading the data
group=group.removeBadStrides; %Removing bad strides
%%

% group=adaptationData.createGroupAdaptData({'ATR01params','ATR02params','ATR03params','ATR04params'});
% group=adaptationData.createGroupAdaptData({'ATS02params','ATS03params','ATS04params','ATS05params','ATS06params','ATS07params',...
%      'ATS08params','ATS09params','ATS10params','ATS11params','ATS12params'});
% group=group.removeBadStrides; %Removing bad strides
age=group.getSubjectAgeAtExperimentDate/12;

%% Define epochs
baseEp=getBaseEpoch;


%Adaptation epochs
% strides=[-150 300 300 300 600];exemptFirst=[0];exemptLast=[0];
strides=[-50 900 300];
% strides=[-40 450 200];

exemptFirst=[1];
exemptLast=[5]; %Strides needed 
names={};
shortNames={};
% cond={'TM Base','Adapt1','Adapt2','Adapt3','Washout'};
% cond={'TM base','gradual adaptation','TM post'}; %Conditions for this group 
cond={'TM base','Adaptation','Post 1'}; %Conditions for this group 
% ep=defineEpochs(cond,cond,strides,exemptFirst,exemptLast,'nanmedian',{'B','A1','A2','A3','P'});

if contains(groupID,'NTS') || contains(groupID,'NTR') || contains(groupID,'CTS') || contains(groupID,'CTR')
    epLong=defineEpochNimbusShoes_longProtocol('nanmean'); 
    baseEp=defineReferenceEpoch('OGNimbus',epLong);
    cond={'TR base','Adaptation','Post 1'}; %Conditions for this group 
    strides=[-50 600 150];
end
ep=defineEpochs(cond,cond,strides,exemptFirst,exemptLast,'nanmean',{'Base','Adapt','Post1'}); %epochs 
%% Define params we care about:
% mOrder={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT', 'VM', 'VL', 'RF', 'HIP','TFL', 'GLU'};
mOrder={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT', 'VM', 'VL', 'RF','TFL', 'GLU'};
% mOrder={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT', 'VM', 'VL', 'RF', 'HIP', 'ADM', 'TFL', 'GLU'};
nMusc=length(mOrder);
type='s';
labelPrefix=fliplr([strcat('f',mOrder) strcat('s',mOrder)]); %To display
labelPrefixLong= strcat(labelPrefix,['_' type]); %Actual names

% %Adding alternative normaliza tion parameters:
l2=group.adaptData{1}.data.getLabelsThatMatch('^Norm');
group=group.renameParams(l2,strcat('N',l2)).normalizeToBaselineEpoch(labelPrefixLong,baseEp,true); %Normalization to max=1 but not min=0
% 
% %Renaming normalized parameters, for convenience:
ll=group.adaptData{1}.data.getLabelsThatMatch('^Norm');
l2=regexprep(regexprep(ll,'^Norm',''),'_s','s');
group=group.renameParams(ll,l2);
newLabelPrefix=strcat(labelPrefix,'s');

%% get data:
padWithNaNFlag=true;
[dataEMG,labels,allDataEMG]=group.getPrefixedEpochData(newLabelPrefix,ep,padWithNaNFlag);
%Flipping EMG:
for i=1:length(allDataEMG)
    aux=reshape(allDataEMG{i},size(allDataEMG{i},1),size(labels,1),size(labels,2),size(allDataEMG{i},3));
    allDataEMG{i}=reshape(flipEMGdata(aux,2,3),size(aux,1),numel(labels),size(aux,4));
end

[~,~,dataContribs]=group.getEpochData(ep,{'netContributionNorm2'},padWithNaNFlag);

%% Save to hdf5 format for sharing with non-Matlab users
EMGdata=cell2mat(allDataEMG);
name='dynamicsData_ATR_NO_HIP.h5';
h5create(name,'/EMGdata',size(EMGdata))
h5write(name,'/EMGdata',EMGdata)
SLA=squeeze(cell2mat(dataContribs));
h5create(name,'/SLA',size(SLA))
h5write(name,'/SLA',SLA)
speedDiff=[zeros(1,50),ones(1,900),zeros(1,300)];
% speedDiff=[zeros(1,40),ones(1,450),zeros(1,200)];
% speedDiff=[zeros(1,40),ones(1,900),zeros(1,200)];
h5create(name,'/speedDiff',size(speedDiff))
h5write(name,'/speedDiff',speedDiff)
% breaks=[zeros(1,150),1,zeros(1,299),1,zeros(1,299),1,zeros(1,299),1,zeros(1,599)];
breaks=[zeros(1,length(speedDiff))];
h5create(name,'/breaks',size(breaks))
h5write(name,'/breaks',breaks)
hdf5write(name,'/labels',l2,'WriteMode','append')
