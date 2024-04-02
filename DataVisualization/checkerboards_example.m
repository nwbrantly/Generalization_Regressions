close all;clear all; clc

groupID ='NTS'; %name of the groups that you want to plot

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

[normalizedGroupData, newLabelPrefix,n_subjects]=creatingGroupdataWnormalizedEMG(groupID,0,refEpForNormalization); % Get normalized data as defult is normalizing by TM base 

%% In this section you define the epoch that you want to plot. 
% Youn may need to costum mode your own defineRegressors function or create
% the epochs as done for the refEpForNormalization

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

EpochsOfInteres={OGbase,TRbase,Pos1_early};


%% Plot inputs 

plotIndSubjects=0; % If you want tp plot individual subjects  
plotGroup=1; %If you want to plot the group average 

%NOTE removebias is 1 you need to define your refence 
removeBias=1; % We are choosing to remove bias 
ref=OGbase; %This is the epoch that we are going to use as reference 
if removeBias==1
    if isempty(ref)
        warning('You need to define your reference epoch')
        return
    end
end
 

if plotIndSubjects
    
    Data= plotEpochsPlusNorm(EpochsOfInteres,normalizedGroupData,newLabelPrefix,plotIndSubjects,0,removeBias,ref);
    
elseif  plotGroup
    
    Data= plotEpochsPlusNorm(EpochsOfInteres,normalizedGroupData,newLabelPrefix,0,plotGroup,removeBias,ref);
    
end

