%% Co-contraction index 
% Script to get the level of co-contraction throghout adaptation. 
% SEe Bank et al, 2017 for reference
clear;clc; close all
%% 1: load and prep data
subID= 'PATR06';
load([subID, 'params.mat'])



%% 2:  EMG normalization of baseline


muscleOrder={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT', 'VM', 'VL', 'RF', 'TFL', 'GLU', 'HIP'};

n_muscles = length(muscleOrder);

ep=defineEpochs_regressionYA('nanmean');
refEp= defineReferenceEpoch('TM base',ep);


newLabelPrefix = defineMuscleList(muscleOrder);

adaptData = adaptData.normalizeToBaselineEpoch(newLabelPrefix,refEp);

ll=adaptData.data.getLabelsThatMatch('^Norm');
l2=regexprep(regexprep(ll,'^Norm',''),'_s','s');
adaptData=adaptData.renameParams(ll,l2);
newLabelPrefix = regexprep(newLabelPrefix,'_s','s');

%% Getting EMG data

% Defining needed variables
data=[];
temp=[];
aux1=[];

Subj = adaptData; %Dummy variable


MusclesofInterest = {'sTAs', 'fTAs','fMGs','sMGs'}; %labels in group ID will be removed for all regression and AE computations;

for i = 1:numel(MusclesofInterest) %loop on the all the muscles
    
    DataIdx=find(contains(Subj.data.labels, {[ MusclesofInterest{i}, ' ']})); %Find data index (row where the muscles are)
    
    if length(DataIdx)<12 % In case the code does not grab all the muscles
        %(It should be 12 gaits phases of the gait cycle)
        DataIdxlast=DataIdx(end)+[1:3];
        DataIdx= [DataIdx; DataIdxlast'];
    end
    
    
    eval([ MusclesofInterest{i},' =  Subj.data.Data(:,DataIdx);']) %Concatenating all the muscle data

end

%% Co-contration index by bin 

CI_sMG = ((2 * sMGs) ./  (sMGs + sTAs))* 100;

CI_sTA = ((2 * sTAs) ./  (sMGs + sTAs))* 100;

CI_fMG = ((2 * fMGs) ./  (fMGs + fTAs))* 100;

CI_fTA = ((2 * fTAs) ./  (fMGs + fTAs))* 100;

