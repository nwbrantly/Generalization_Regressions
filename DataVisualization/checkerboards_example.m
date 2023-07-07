groupID ='BATR'; %name of the groups that you want to plot 

[normalizedGroupData, newLabelPrefix,n_subjects]=creatingGroupdataWnormalizedEMG(groupID); % Get normalized data as defult is normalizing by TM base 

ep=defineRegressorsDynamicsFeedback('nanmean'); %Defining epochs of interest. This needs to be costum to your study 

% Epochs that you want to plot 
OGbase=defineReferenceEpoch('OG base',ep);
TMbase= defineReferenceEpoch('TM base',ep);
Pos1_Late=defineReferenceEpoch('Post1_{Early}',ep);
Pos2_Late=defineReferenceEpoch('Post2_{Early}',ep);

EpochsOfInteres={OGbase,TMbase,Pos1_Late,Pos2_Late};


plotIndSubjects=0; % If you want tp plot individual subjects  
plotGroup=1; %If you want to plot the group average 

%NOTE removebias is 1 you need to define your refence 

removeBias=1; % We are choosing to remove bias 
ref=[]; %This is the epoch that we are going to use as reference 

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
    plotEpochsPlusNorm(EpochsOfInteres,normalizedGroupData,newLabelPrefix,0,plotGroup,removeBias,ref)
    
end

