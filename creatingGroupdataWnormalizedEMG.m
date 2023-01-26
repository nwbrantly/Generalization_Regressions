function [normalizedGroupData, newLabelPrefix]=creatingGroupdataWnormalizedEMG(groupID)
%This function is to created the group data and also normalize the EMG. In
%this version I am normalizing to TM base (tied belt right before
%adpatation)

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
n_muscles = length(muscleOrder);

ep=defineRegressorsDynamicsFeedback('nanmean');
refEpTM = defineReferenceEpoch('TM base',ep);


GroupData=adaptationData.createGroupAdaptData(sub); %loading the data
GroupData=GroupData.removeBadStrides; %Removing bad strides

newLabelPrefix = defineMuscleList(muscleOrder);
normalizedGroupData = GroupData.normalizeToBaselineEpoch(newLabelPrefix,refEpTM); 
ll=normalizedGroupData.adaptData{1}.data.getLabelsThatMatch('^Norm');
l2=regexprep(regexprep(ll,'^Norm',''),'_s','s');
normalizedGroupData=normalizedGroupData.renameParams(ll,l2);
newLabelPrefix = regexprep(newLabelPrefix,'_s','s');

end