%%
clear; close all; clc;

%%

%% Defining the epoch that we are going to use for data normalization
epochNames={'TM base'};
condition= {'TM base'}; %Change conditions names to your own! 
strideNo=[-40]; %Positive vaues define inital; negative values define # strides at end of that condition
exemptFirst=0; %Number of strides you want to ignore at the beginning of the condition
exemptLast=5; %Number of strides you want to ignore at the end of the condition
summaryMethod={'nanmean'}; %Method to analyze bar plots 
shortName=[];
[refEpForNormalization] = defineEpochs(epochNames,condition,strideNo,exemptFirst,exemptLast,summaryMethod,shortName);

%% 
groupID ='C3S01'; %Group of interest 
[group, newLabelPrefix,n,subID]=creatingGroupdataWnormalizedEMG(groupID,1,refEpForNormalization); % Creating the groupData normalized
sesion1=1;
adaptation=0;
negative=0;
if strcmp(groupID,'C3')
    sesion1  = questdlg('Is this session 1?', ...
        'Session of interest', ...
        'Yes','No','None');
    
    if strcmpi(sesion1, 'Yes')
        sesion1 = 1;
    else
        sesion1 = 2;
    end
end

%% Removing bad muscles 
%This script make sure that we always remove the same muscle for the
%different analysis 
removeBadmuscles=1;
if removeBadmuscles==1
group= RemovingBadMuscleToSubj(group);
end
%% Define epochs depending on the group data

if contains(groupID,'BAT')

    if adaptation==1
        strides=[-40 450]; %Number per strides per condition
        cond={'TM base','Adaptation'}; %Conditions for this group
    elseif  contains(groupID,'BATR')
        strides=[-40 300]; %Number per strides per condition
        cond={'TM base','Post 1'}; %Conditions for this group
    elseif  contains(groupID,'BATS')
        strides=[-40 300]; %Number per strides per condition
        cond={'OG base','Post 1'}; %Conditions for this group
        
    end
    exemptFirst=[1]; %ignore inital strides
    exemptLast=[5]; %Strides needed
    

elseif contains(groupID,'CTS') || contains(groupID,'NTS') || contains(groupID,'VATS')
    strides=[-40 150]; %Number per strides per condition
    cond={'OG base','Post 1'}; %Conditions for this group
    exemptFirst=[1]; %ignore inital strides
    exemptLast=[5]; %Strides needed
    
elseif contains(groupID,'CTR') || contains(groupID,'NTR') || contains(groupID,'VATR')
     strides=[-40 150]; %Number per strides per condition
    cond={'TR base','Post 1'}; %Conditions for this group
    exemptFirst=[1]; %ignore inital strides
    exemptLast=[5]; %Strides needed

elseif contains(groupID,'C3')
    
    if  adaptation==1 
        strides=[-40 900]; %Number per strides per condition
        cond={'TM base','Adaptation'}; %Conditions for this group
    else
    if sesion1==1    
        cond={'OG base','Post 1'}; %Conditions for this group
    elseif sesion1==2
        cond={'TM base','Post 1'}; %Conditions for this group
    end
    
    strides=[-40 350]; %Number per strides per condition

    end

    exemptFirst=[1]; %ignore inital strides
    exemptLast=[5]; %Strides needed   
    
elseif contains(groupID,'MWS')
    if  adaptation==1
        strides=[-40 900]; %Number per strides per condition
        cond={'NIM base','Adaptation'}; %Conditions for this group
    elseif negative==1
        strides=[-40 900]; %Number per strides per condition
          cond={'NIM base','OG Tied'};  %Conditions for this group
           strides=[-40 20]; %Number per strides per condition
    
    else
        cond={'OG base','Post 1'}; %Conditions for this group
        
        strides=[-40 200]; %Number per strides per condition
        
    end
    
    exemptFirst=[1]; %ignore inital strides
    exemptLast=[5]; %Strides needed
    
end

ep=defineEpochs(cond,cond,strides,exemptFirst,exemptLast,'nanmean',{'Base','Post1'}); %epochs
%% Pick muscles that you wanted to get the data from 
%%%% load and prep data
% muscleOrder={'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT'};
muscleOrder={'TA','PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT', 'VM', 'VL', 'RF', 'HIP','TFL', 'GLU'};
newLabelPrefix2 = defineMuscleList(muscleOrder); %List of muscle
newLabelPrefix2 = regexprep(newLabelPrefix2,'_s','s'); %removing the underscore "_"

for m=1:length(newLabelPrefix2)

    idx(m)=find(strcmp(newLabelPrefix,newLabelPrefix2(m))); %getting the index of the muscles

end

 wanted_Muscles= newLabelPrefix(sort(idx)); %It needs to be the muscle form the slow leg first
 newLabelPrefix= wanted_Muscles;
 
%% get data:
padWithNaNFlag=true;
[dataEMG,labels,allDataEMG2]=group.getPrefixedEpochData(newLabelPrefix,ep,padWithNaNFlag); 

%Flipping EMG:
for i=1:length(allDataEMG2)
    aux=reshape(allDataEMG2{i},size(allDataEMG2{i},1),size(labels,1),size(labels,2),size(allDataEMG2{i},3));
    allDataEMG2{i}=reshape(flipEMGdata(aux,2,3),size(aux,1),numel(labels),size(aux,4));
end

[~,~,dataContribs]=group.getEpochData(ep,{'netContributionNorm2'},padWithNaNFlag);


%% Getting the regressors values
%% 
if contains(groupID,'BAT')
    ep=defineRegressorsDynamicsFeedback('nanmean');
    epochOfInterest={'Ramp','PosShort_{late}','Adaptation_{early}','Adaptation','Optimal','NegShort_{late}','NegShort_{early}','TM base'};
elseif contains(groupID,'CTS') || contains(groupID,'CTR') || contains(groupID,'NTS') || contains(groupID,'NTR') || contains(groupID,'VATS') || contains(groupID,'VATR')
    ep=defineEpochNimbusShoes('nanmean');
    epochOfInterest={'SplitNeg','Adaptation','SplitPos','TRbase','OGbase'};
elseif contains(groupID,'C3') 
    ep=defineRegressors_StrokeC3('nanmean');
    epochOfInterest={'PosShort_{late}','Adaptation_{early}','Adaptation','NegShort_{late}','TM base','OG base'};
elseif contains(groupID,'MWS') 
   ep= defineEpochMW('nanmean');
   epochOfInterest={'NIMbase','OGbase','Adaptation_{early}','Adaptation_{late}','NegShort','Post1_{early}'}; 
end


    padWithNaNFlag=true; %If no enough steps fill with nan, let this on
for l=1:length(epochOfInterest)
   
    ep2=defineReferenceEpoch(epochOfInterest{l},ep);
    
    [dataEMG,labels,allDataEMG]=group.getPrefixedEpochData(newLabelPrefix,ep2,padWithNaNFlag); %Getting the data
    
    %Flipping EMG:0
    for i=1:length(allDataEMG)
        aux=reshape(allDataEMG{i},size(allDataEMG{i},1),size(labels,1),size(labels,2),size(allDataEMG{i},3));
        allDataEMG{i}=reshape(flipEMGdata(aux,2,3),size(aux,1),numel(labels),size(aux,4));
        
    end
    
  regressors{l}=nanmean(allDataEMG{:},1);


end

%% Reorganize data
EMGdata=cell2mat(allDataEMG2);
regressors=cell2mat(regressors');

%% Save to hdf5 format for sharing with non-Matlab users

if adaptation==1
    name=['dynamicsData_',groupID,'_subj_', num2str(n),'_Session_',num2str(sesion1),'_RemoveBadMuscles_', num2str(removeBadmuscles),'_',datestr(now,'dd-mmmm-yyyy'),'_Adaptation','.h5'];
else
    name=['dynamicsData_',groupID,'_subj_', num2str(n),'_Session_',num2str(sesion1),'_RemoveBadMuscles_', num2str(removeBadmuscles),'_',datestr(now,'dd-mmmm-yyyy'),'_Post-Adaptation','.h5'];
    
end

h5create(name,'/EMGdata',size(EMGdata))
h5write(name,'/EMGdata',EMGdata)

h5create(name,'/Regressors',size(regressors))
h5write(name,'/Regressors',regressors)

hdf5write(name,'/SubID',subID(:),'WriteMode','append')

hdf5write(name,'/Epochs',epochOfInterest(:),'WriteMode','append')

SLA=squeeze(cell2mat(dataContribs));
h5create(name,'/SLA',size(SLA))
h5write(name,'/SLA',SLA)
% speedDiff=[zeros(1,abs(strides(1))),ones(1,strides(2)),zeros(1,(strides(3)))];
speedDiff=[zeros(1,abs(strides(1))),ones(1,strides(2))];
h5create(name,'/speedDiff',size(speedDiff))
h5write(name,'/speedDiff',speedDiff)
breaks=[zeros(1,length(speedDiff))];
h5create(name,'/breaks',size(breaks))
h5write(name,'/breaks',breaks)
hdf5write(name,'/labels',newLabelPrefix(:),'WriteMode','append')
