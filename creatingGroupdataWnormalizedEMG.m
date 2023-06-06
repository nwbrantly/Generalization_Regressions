function [GroupData, newLabelPrefix,n_subjects]=creatingGroupdataWnormalizedEMG(groupID,invMuscles)
%This function is to created the group data and also normalize the EMG. In
%this version I am normalizing to TM base (tied belt right before
%adpatation). We are normalizing everthing by TM base for consistance

files = dir ([groupID '*params.mat']);


n_subjects = size(files,1);

ii=0;
for i =1:n_subjects
    ii=1+ii;
    sub{ii} = files(i).name;
    subID{ii} = sub{ii}(1:end-10);
end

subID % This display what are the ID of the partocopants that you are using. 

%%

%%%% load and prep data
muscleOrder={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT', 'VM', 'VL', 'RF', 'HIP','TFL', 'GLU'};
%%%% load and prep data Posterior muscles
% muscleOrder={ 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT'};
newLabelPrefix = defineMuscleList(muscleOrder); %List of muscle 
n_muscles=length(muscleOrder)*2;


if nargin>1 
    if invMuscles==1 %For the prepocess code this needs to run bc we need to organize the data correctly 
    newLabelPrefix=newLabelPrefix([n_muscles:-1:n_muscles/2+1 n_muscles/2:-1:1]);
    end
end

ep=defineRegressorsDynamicsFeedback('nanmedian'); %loading the variables for multiple epochs of interest
refEpTM = defineReferenceEpoch('TM base',ep); %extracting the info for the epoch that want to use to normalize the data
 
GroupData=adaptationData.createGroupAdaptData(sub); %loading the data
GroupData=GroupData.removeBadStrides; %Removing bad strides


GroupData = GroupData.normalizeToBaselineEpoch(newLabelPrefix,refEpTM,true); %Normalizing the data
ll=GroupData.adaptData{1}.data.getLabelsThatMatch('^Norm'); %Getting the label witht the word Norm 
l2=regexprep(regexprep(ll,'^Norm',''),'_s','s'); %Removing the word norm. This is just to make our life easier
GroupData=GroupData.renameParams(ll,l2); %Rename the variable to not have the word Norm 
newLabelPrefix = regexprep(newLabelPrefix,'_s','s'); %removing the underscore "_"

end