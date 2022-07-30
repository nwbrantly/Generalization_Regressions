%% Co-contraction index 
% Script to get the level of co-contraction throghout adaptation. 
% See Bank et al, 2017 for reference

clear; close all; clc;

% set script parameters, SHOULD CHANGE/CHECK THIS EVERY TIME.
groupID = 'BATR';
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
    
    strides=[-40 450];
    cond={'TM base', 'Adaptation'}; %Conditions for this group
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
ep=defineEpochs(cond,cond,strides,exemptFirst,exemptLast,'nanmean',{'TM base','Adapt'}); %epochs 
%% Define params we care about:
% mOrder={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT', 'VM', 'VL', 'RF', 'HIP','TFL', 'GLU'};
mOrder={'TA', 'MG'};
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
EMGdata=cell2mat(allDataEMG);
muscPhaseIdx=1:size(EMGdata,2);%336; %14 muscles per leg %360; %All muscles
Y=EMGdata(:,muscPhaseIdx,:);
Y=nanmedian(Y,3); %Median across subjs
idx=1:12:49;
sMG=Y(:,idx(1):idx(2)-1);
sTA=Y(:,idx(2):idx(3)-1);
fMG= Y(:,idx(3):idx(4)-1);
fTA= Y(:,idx(4):idx(5)-1);


%% Co-contration index by bin 

%﻿During LR, ISw, MSw, and TSw, the TA should be the agonist muscle, 
%while the MG should be the antagonist: CI= I_ant/(I_ant+I_ago)
% MG antagonsitic for 2,3,4,5,6,7
CI_fMG = ((2 * fMG) ./  (fMG + fTA))* 100;
CI_sMG = ((2 * sMG) ./  (sMG + sTA))* 100;



% ﻿During MSt and TSt, the MG should be the agonist and the TA the
% antagonist muscle:
%TA antagonistic for 1,9,10,11,12
CI_sTA = ((2 * sTA) ./  (sMG + sTA))* 100;
CI_fTA = ((2 * fTA) ./  (fMG + fTA))* 100;

%%
figure
plot(movmean(CI_sMG(:,4),5))
hold on 
plot(movmean(CI_fMG(:,4),5))
title('Co-contraction index = 100 * (2 * MG)/ (MG + TA)')
legend('slow','fast')
ylabel('CI during mid stance (subin=4) ')
set(gcf,'color','w')
xlabel('strides')