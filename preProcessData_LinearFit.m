%%
%Running this script requires labTools (github.com/pittSMLlab/labTools/)
% addpath(genpath('../../../splitbelt-EMG-adaptation/'))
% addpath(genpath('../../../EMG-LTI-SSM/'))
clear; close all; clc;

% set script parameters, SHOULD CHANGE/CHECK THIS EVERY TIME.
groupID = 'BATS';
scriptDir = cd;% fileparts(matlab.desktop.editor.getActiveFilename);
files = dir ([groupID '*params.mat']);



sub = {};
subID = {};


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
age=group.getSubjectAgeAtExperimentDate/12;

%% Define epochs
baseEp=getBaseEpoch;


%Adaptation epochs
% strides=[-150 300 300 300 600];exemptFirst=[0];exemptLast=[0];
% strides=[-50 900 300];
splits=0;
if splits==1
    strides=[-40 980 100];
    cond={'TM base','Multiple pos short splits','TM mid 2'}; %Conditions for this group
else
    
    strides=[-40 440 200];
    cond={'TM base','Adaptation','Post 1'}; %Conditions for this group
%     strides=[50];
%     cond={'Post 1'}; %Conditions for this group
end
exemptFirst=[1];
exemptLast=[5]; %Strides needed 
names={};
shortNames={};


if contains(groupID,'NTS') || contains(groupID,'NTR') || contains(groupID,'CTS') || contains(groupID,'CTR')
    epLong=defineEpochNimbusShoes_longProtocol('nanmean'); 
    baseEp=defineReferenceEpoch('OGNimbus',epLong);
    cond={'TR base','Adaptation','Post 1'}; %Conditions for this group 
    strides=[-50 600 150];
end
ep=defineEpochs(cond,cond,strides,exemptFirst,exemptLast,'nanmean',{'Base','Adapt','Post1'}); %epochs 
% ep=defineEpochs(cond,cond,strides,exemptFirst,exemptLast,'nanmean',{'Post1'}); %epochs 
%% Define params we care about:
mOrder={'TA','PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT', 'VM', 'VL', 'RF', 'HIP','TFL', 'GLU'};
nMusc=length(mOrder);
type='s';
labelPrefix=fliplr([strcat('f',mOrder) strcat('s',mOrder)]); %To display
labelPrefixLong= strcat(labelPrefix,['_' type]); %Actual names

% %Adding alternative normaliza tion parameters:
l2=group.adaptData{1}.data.getLabelsThatMatch('^Norm');
group=group.renameParams(l2,strcat('N',l2)).normalizeToBaselineEpoch(labelPrefixLong,baseEp,true); %Normalization to max=1 but not min=0


% %Renaming normalized parameters, for convenience:
ll=group.adaptData{1}.data.getLabelsThatMatch('^Norm');
l2=regexprep(regexprep(ll,'^Norm',''),'_s','s');
group=group.renameParams(ll,l2);
newLabelPrefix=strcat(labelPrefix,'s');

%% Set bad muscles to nan
 
 removeBadmuscles=0;
if removeBadmuscles==1
    [RemovedData]=RemoveBadMuscles(group,groupID);
    group=RemovedData;
end

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
name=['dynamicsData_',groupID,'_subj_', num2str(size(subID,2)),'_RemoveBadMuscles', num2str(removeBadmuscles),'_splits_',num2str(splits),'.h5'];
h5create(name,'/EMGdata',size(EMGdata))
h5write(name,'/EMGdata',EMGdata)
SLA=squeeze(cell2mat(dataContribs));
h5create(name,'/SLA',size(SLA))
h5write(name,'/SLA',SLA)
speedDiff=[zeros(1,abs(strides(1))),ones(1,strides(2)),zeros(1,(strides(3)))];
% speedDiff=[ones(1,strides(1))];
h5create(name,'/speedDiff',size(speedDiff))
h5write(name,'/speedDiff',speedDiff)
breaks=[zeros(1,length(speedDiff))];
h5create(name,'/breaks',size(breaks))
h5write(name,'/breaks',breaks)
hdf5write(name,'/labels',l2,'WriteMode','append')
