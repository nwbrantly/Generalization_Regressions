%% Adding Norm to groupAdaptationData

% Load data and Find norms for the entire time courses
%This code will find Euclidean norms for the entire time courses
%Created by DMM0 5/2022 - Update by DMMO Oct 2023

% 1) load subjects
% 2) EMG normalization of baseline
% 3) Remove bad muscles making then zero. We are computing the norm
% 4) Computing Stride by stride norm
% 5) Compute bias removed stride by stride norm
% 6) Saving params file 
clear; clc; close all
%% Define conditions
saveData=1 %if you want to save your data 

%% load subjects  and prep data (Step 1 and 2)


epochNames={'TM Base'};
condition= {'TM Base'}; %Change conditions names to your own! 
strideNo=[-40]; %Positive vaues define inital; negative values define # strides at end of that condition
exemptFirst=0; %Number of strides you want to ignore at the beginning of the condition
exemptLast=5; %Number of strides you want to ignore at the end of the condition
summaryMethod={'nanmean'}; %Method to analyze bar plots 
shortName=[];
[refEpForNormalization] = defineEpochs(epochNames,condition,strideNo,exemptFirst,exemptLast,summaryMethod,shortName);

groupID = {'SAH12'};
scriptDir = cd;
[normalizedGroupData, newLabelPrefix,n_subjects,subID]=creatingGroupdataWnormalizedEMG(groupID{1},0,refEpForNormalization);
removeMuscles=1;
%% Removing bad muscles (Step 3) 

%This script make sure that we always remove the same muscle for the
%different analysis. 
% NOTE: THIS IS HARD-CODED PLEASE UPDATE YOUR CODE

if removeMuscles==1
    normalizedGroupData= RemovingBadMuscleToSubj(normalizedGroupData);
end

%% Compute bias, unbias and individual muscle stride by stride norm (Step 4 and 5)
newLabelPrefix = newLabelPrefix(8); % SAH_Pilot02: 8 = sMG, others: 12 = sMG
[normalizedGroupData]=AddingNorm(normalizedGroupData,subID,newLabelPrefix,groupID);

%% SAVE GROUP DATA (Step 6)

    group= normalizedGroupData;

if saveData==1
    save([groupID{1}, '_NormEMG.mat'],'group')
end

%%  This section is to plot the results

group2=[];
group2{1}=group;


conditions = group.adaptData{1}.metaData.conditionName; % {'OG base','TM base','Adaptation',...
    % 'Post 1','Post 2'};


params={'NormEMG','sMGsNorm'}; % 'netContributionNorm2',
poster_colors;
     
binwidth=5; %Window of the running average
trialMarkerFlag=0; %1 if you want to separete the time course by trial 0 to separece by condition 
indivFlag=0; %0 to plot group mean 1 to plot indiv subjects
indivSubs=[]; %Use when you want to plot a specidfic subject in a group 
colorOrder=[];%[p_red; p_orange; p_plum;p_fade_green]; %Let the function take care of this at least you wanted in a specific set of color then by my guess and add the list here
biofeedback= 0; % At least that you are providing with biofeedback to the subject
removeBiasFlag=0; %if you want to remove bias 
%%Groups names 
labels=[];
filterFlag=[]; 
% figure 
% p=subplot(1,1,1);
plotHandles=[];
alignEnd=0; % # strides align at the end of the trial (PLAY with it as see what happens)
alignIni=0; %  # strides align at the beginning of the trial (PLAY with it as see what happens) 

adaptData=cellfun(@(x) x.adaptData,group2,'UniformOutput',false); %Notice that adaptDataGroups(1) decide that I only want to plot the CG group 
[figh,avg,indv]=adaptationData.plotAvgTimeCourse(adaptData,params,conditions,binwidth,trialMarkerFlag,...
    indivFlag,indivSubs,colorOrder,biofeedback,removeBiasFlag,labels,filterFlag,plotHandles,alignEnd,alignIni);
% legend('AutoUpdate','off')
% yline(0)

% yline(nanmean(avg.NormEMGasym.TMbase.trial1(1,end-15:end)))