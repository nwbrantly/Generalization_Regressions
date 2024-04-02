close all;clear all; clc

groupID ='NTS'; %name of the groups that you want to plot 

[normalizedGroupData, newLabelPrefix,n_subjects]=creatingGroupdataWnormalizedEMG(groupID); % Get normalized data as defult is normalizing by TM base 

% ep=defineRegressors_StrokeC3('nanmean'); %Defining epochs of interest. This needs to be costum to your study 

ep=defineEpochNimbusShoes('nanmean'); 
% Epochs that you want to plot 
OGbase=defineReferenceEpoch('OGbase',ep);
TRbase= defineReferenceEpoch('TRbase',ep);
% EarlyAdapt= defineReferenceEpoch('Adaptation_{early}',ep);
LateAdapt= defineReferenceEpoch('Adaptation',ep);
Pos1_early=defineReferenceEpoch('Post1_{Early}',ep);
% Pos2_early=defineReferenceEpoch('Post2_{Early}',ep);
% % Optimal=defineReferenceEpoch('Optimal',ep);
NegShort=defineReferenceEpoch('SplitNeg',ep);
% PostShort= defineReferenceEpoch('PosShort_{late}',ep);
% EpochsOfInteres={OGbase,TMbase,EarlyAdapt,LateAdapt,Pos1_Late,Pos2_Late};
%%
% EpochsOfInteres={TRbase,NegShort,LateAdapt};
EpochsOfInteres={OGbase,TRbase,Pos1_early};
plotIndSubjects=0; % If you want tp plot individual subjects  
plotGroup=1; %If you want to plot the group average 

%NOTE removebias is 1 you need to define your refence 

removeBias=1; % We are choosing to remove bias 
ref=OGbase; %This is the epoch that we are going to use as reference 
%  ref=TRbase; %This is the epoch that we are going to use as reference 
if removeBias==1
    if isempty(ref)
        warning('You need to define your refence')
        return
    end
end

if plotIndSubjects
    plotEpochsPlusNorm(EpochsOfInteres,normalizedGroupData,newLabelPrefix,plotIndSubjects,0,removeBias,ref)
end

if plotGroup
    Data= plotEpochsPlusNorm(EpochsOfInteres,normalizedGroupData,newLabelPrefix,0,plotGroup,removeBias,ref)
    
end

