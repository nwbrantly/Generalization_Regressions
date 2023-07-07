function [Data,regressorNames,fh]=getRegressionsData(Labels,session1,session2,plotIndSubjects, plotGroup,NegShort,TMbeforeNeg,PostShort,TMbeforePost,AdaptLate,Post1Early,Post1Late,Post2Early, OGpostPosEarly, OGbase, EnvBase,subIdx,flip)

if nargin<17 || isempty(subIdx)
    subIdx=[];
end

n_subjects=length(session1.ID);
regModelVersion = 'default';
regressorNames = {'Adapt','Noadapt','MultiContextSwitch','Trans1','Trans2','Trans3'};

if  plotIndSubjects
    
        disp(['subject=', num2str(subIdx)])
        
        
        adaptDataSubject = session1.adaptData{1, subIdx};
        
        if ~isempty(session2)
            adaptDataSubjectSession2 = session2.adaptData{1, subIdx};
        else
            adaptDataSubjectSession2 = session1.adaptData{1, subIdx};
        end
          
        
        Data = {}; %in order: adapt, dataEnvSwitch, dataTaskSwitch, dataTrans1, dataTrans2
        %all labels should be the same, no need to save again.

        
        
        regressorNames = {'Adapt','Noadapt','MultiContextSwitch','Trans1','Trans2','Trans3'};
        
        [Data{1}]=adaptDataSubjectSession2.getCheckerboardsData(Labels,NegShort,TMbeforeNeg,flip);

        [Data{2}] = adaptDataSubjectSession2.getCheckerboardsData(Labels,TMbeforePost,PostShort,flip);

        [Data{3}] = adaptDataSubject.getCheckerboardsData(Labels,OGbase,EnvBase,flip);
          
        [Data{4}] = adaptDataSubject.getCheckerboardsData(Labels,Post1Early,AdaptLate,flip); %Post1 early - Adaptation_{SS} , transition 1

        [Data{5}] = adaptDataSubject.getCheckerboardsData(Labels,Post2Early,Post1Late,flip); %Post2 early - Post 1 late, transition 2

        [Data{6}] = adaptDataSubjectSession2.getCheckerboardsData(Labels,OGpostPosEarly,PostShort,flip); % OGpost Pos - Pos Short


end

if  plotGroup
    
    regModelVersion = 'default';

    if ~isempty(session2)
        session2Data = session2;
    else
        session2Data = session1;
    end
   
       summaryflag='nanmedian';     
        
        Data = {}; %in order: {'Adapt','WithinContextSwitch','MultiContextSwitch','Trans1','Trans2'};
      
        [Data{1}]=session2Data.getCheckerboardsData(Labels,NegShort,TMbeforeNeg,flip,summaryflag); %  EMG_split(-) - TM base slow, adaptation

        [Data{2}] = session2Data.getCheckerboardsData(Labels,TMbeforePost,PostShort,flip,summaryflag); % Noadapt (env-driven/within-env), - EMGon(+) = TM base - EMG_on(+)
    
        [Data{3}] = session1.getCheckerboardsData(Labels,OGbase,EnvBase,flip,summaryflag); %  OG base - TR base = -(TR base - OG base), env switching
 
        [Data{4}] = session1.getCheckerboardsData(Labels,Post1Early,AdaptLate,flip,summaryflag);%Post1 early - Adaptation_{SS}, transition 1
        
        [Data{5}] = session1.getCheckerboardsData(Labels,Post2Early,Post1Late,flip,summaryflag); %TM post VR early - OG post late, transition 2
        
        [Data{6}] = session2Data.getCheckerboardsData(Labels,OGpostPosEarly,PostShort,flip,summaryflag); %TM post VR early - OG post late, transition 2

        
    
end


end