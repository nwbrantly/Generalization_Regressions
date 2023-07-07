groupID ='BATR';
[normalizedGroupData, newLabelPrefix,n_subjects]=creatingGroupdataWnormalizedEMG(groupID);

ep=defineRegressorsDynamicsFeedback('nanmean');

OGbase=defineReferenceEpoch('OG base',ep);
TMbase= defineReferenceEpoch('TM base',ep);
Pos1_Late=defineReferenceEpoch('Post1_{early}',ep);
Pos2_Late=defineReferenceEpoch('Post2_{early}',ep);

EpochsOfInteres={OGbase,TMbase,Pos1_Late,Pos2_Late};

plotIndSubjects=1
 plotGroup=1

if plotIndSubjects
    plotEpochsPlusNorm(EpochsOfInteres,normalizedGroupData,newLabelPrefix,1,0)
end

if plotGroup
    plotEpochsPlusNorm(EpochsOfInteres,normalizedGroupData,newLabelPrefix,0,1)
    
end

