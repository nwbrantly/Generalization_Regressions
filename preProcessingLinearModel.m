%%
clear; close all; clc;
%%
groupID ='BATS'; %Group of interest 
[group2, newLabelPrefix,n]=creatingGroupdataWnormalizedEMG(groupID,1); % Creating the gropData normalized
%% Removing bad muscles 
%This script make sure that we always remove the same muscle for the
%different analysis 
removeBadmuscles=0;
if removeBadmuscles==1
group= RemovingBadMuscleToSubj(group);
end
%% Define epochs
splits=0;
baseEp=getBaseEpoch;
strides=[-40 440 200];
cond={'TM base','Adaptation','Post 1'}; %Conditions for this group

exemptFirst=[1];
exemptLast=[5]; %Strides needed
names={};
shortNames={};


ep=defineEpochs(cond,cond,strides,exemptFirst,exemptLast,'nanmean',{'Base','Adapt','Post1'}); %epochs
%% get data:
padWithNaNFlag=true;
[dataEMG,labels,allDataEMG2]=group2.getPrefixedEpochData(newLabelPrefix,ep,padWithNaNFlag);

%Flipping EMG:
for i=1:length(allDataEMG2)
    aux=reshape(allDataEMG2{i},size(allDataEMG2{i},1),size(labels,1),size(labels,2),size(allDataEMG2{i},3));
    allDataEMG2{i}=reshape(flipEMGdata(aux,2,3),size(aux,1),numel(labels),size(aux,4));
end

[~,~,dataContribs]=group2.getEpochData(ep,{'netContributionNorm2'},padWithNaNFlag);

%% Save to hdf5 format for sharing with non-Matlab users
EMGdata=cell2mat(allDataEMG2);
name=['dynamicsData_',groupID,'_subj_', num2str(n),'_RemoveBadMuscles', num2str(removeBadmuscles),'_splits_',num2str(splits),'V3','.h5'];
h5create(name,'/EMGdata',size(EMGdata))
h5write(name,'/EMGdata',EMGdata)
SLA=squeeze(cell2mat(dataContribs));
h5create(name,'/SLA',size(SLA))
h5write(name,'/SLA',SLA)
speedDiff=[zeros(1,abs(strides(1))),ones(1,strides(2)),zeros(1,(strides(3)))];
h5create(name,'/speedDiff',size(speedDiff))
h5write(name,'/speedDiff',speedDiff)
breaks=[zeros(1,length(speedDiff))];
h5create(name,'/breaks',size(breaks))
h5write(name,'/breaks',breaks)
hdf5write(name,'/labels',newLabelPrefix(:),'WriteMode','append')
