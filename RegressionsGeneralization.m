function [Data,regressorNames,fh]=RegressionsGeneralization(Labels,session1,session2,plotIndSubjects, plotGroup,NegShort,TMbeforeNeg,PostShort,TMbeforePost,AdaptLate,Post1Early,Post1Late,Post2Early, OGpostPosEarly, OGbase, EnvBase,subIdx,flip)

if nargin<17 || isempty(subIdx)
    subIdx=[];
end

n_subjects=length(session1.ID);

if  plotIndSubjects
    
        disp(['subject=', num2str(subIdx)])
        
        
        adaptDataSubject = session1.adaptData{1, subIdx};
        
        if ~isempty(session2)
            adaptDataSubjectSession2 = session2.adaptData{1, subIdx};
        else
            adaptDataSubjectSession2 = session1.adaptData{1, subIdx};
        end
        
        fh=figure('Units','Normalized','OuterPosition',[0 0 1 1],'NumberTitle', 'off', 'Name',session1.ID{subIdx});
        ph=tight_subplot(1,6,[.03 .005],.04,.04);
        
        
        
        Data = {}; %in order: adapt, dataEnvSwitch, dataTaskSwitch, dataTrans1, dataTrans2
        %all labels should be the same, no need to save again.
        
        
        regressorNames = {'MultiContextAdapt','EnvTransition','MultiContextSwitch','Trans1','Trans2','Trans3'};
        
        [~,~,labels,Data{1},dataRef2]=adaptDataSubjectSession2.plotCheckerboards(Labels,NegShort,fh,ph(1,1),TMbeforeNeg,flip); %  EMG_split(-) - TM base slow, adaptation
        title('Within-Env-Adapt: NegShort-TMSlow') % (-) early - TM slow
        regressorNames{1} = 'Adapt';
        
        [~,~,~,Data{2},~] = adaptDataSubjectSession2.plotCheckerboards(Labels,TMbeforePost,fh,ph(1,2),PostShort,flip); % Noadapt (env-driven/within-env), - EMGon(+) = TM base - EMG_on(+)
        title('Multi-Env-Switch: TMfast-PosShort')
        regressorNames{2} = 'Noadapt';
        
        [~,~,~,Data{3},~] = adaptDataSubject.plotCheckerboards(Labels,OGbase,fh,ph(1,3),EnvBase,flip);
        title('Env-Switch: OG-OGNimbus')  % OG - OGNimbus
        regModelVersion = 'default'
        
        [~,~,~,Data{4},~] = adaptDataSubject.plotCheckerboards(Labels,Post1Early,fh,ph(1,4),AdaptLate,flip); %Post1 early - Adaptation_{SS} , transition 1
        title('Transition 1')
        [~,~,~,Data{5},~] = adaptDataSubject.plotCheckerboards(Labels,Post2Early,fh,ph(1,5),Post1Late,flip); %Post2 early - Post 1 late, transition 2
        title('Transition 2')
        
        [~,~,~,Data{6},~] = adaptDataSubjectSession2.plotCheckerboards(Labels,OGpostPosEarly,fh,ph(1,6),PostShort,flip); % OGpost Pos - Pos Short
        title('Transition 3: Short Split')
        
        
        
        
        set(ph(:,1),'CLim',[-1 1]*1);
        set(ph(:,2:end),'YTickLabels',{},'CLim',[-1 1]*1);
        set(ph,'FontSize',8)
        pos=get(ph(1,end),'Position');
        axes(ph(1,end))
        colorbar
        set(ph(1,end),'Position',pos);
        set(gcf,'color','w');
        

end

if  plotGroup

    if ~isempty(session2)
        session2Data = session2;
    else
        session2Data = session1;
    end
    
%     flip = [1];
       summaryflag='nanmedian';

        fh=figure('Units','Normalized','OuterPosition',[0 0 1 1]);
        ph=tight_subplot(1,6,[.03 .005],.04,.04);
        
        
        Data = {}; %in order: {'Adapt','WithinContextSwitch','MultiContextSwitch','Trans1','Trans2'};
        %                         [~,~,labels,dataE{1},dataRef{1}]=session1.plotCheckerboards(Labels,ep,fh,ph(1,2:end),refEp,flip);%Second, the rest:
       
        
        [~,~,labels,Data{1},dataRef2]=session2Data.plotCheckerboards(Labels,NegShort,fh,ph(1,1),TMbeforeNeg,flip,summaryflag); %  EMG_split(-) - TM base slow, adaptation
        d = nanmedian(Data{1}, 4);
        title(['Within-Env-Adapt: NegShort-TMSlow'] ,[ 'Norm=', num2str(norm(reshape(d,[],1)))]) % (-) early - TM slow
        regressorNames{1} = 'Adapt';
        
        [~,~,~,Data{2},~] = session2Data.plotCheckerboards(Labels,TMbeforePost,fh,ph(1,2),PostShort,flip,summaryflag); % Noadapt (env-driven/within-env), - EMGon(+) = TM base - EMG_on(+)
        d = nanmedian(Data{2}, 4);
        title(['Multi-Env-Switch: TMfast-PosShort'] ,[ 'Norm=', num2str(norm(reshape(d,[],1)))])
        regressorNames{2} = 'WithinContextSwitch';
        
        [~,~,~,Data{3},~] = session1.plotCheckerboards(Labels,OGbase,fh,ph(1,3),EnvBase,flip,summaryflag); %  OG base - TR base = -(TR base - OG base), env switching
        d = nanmedian(Data{3}, 4);
        title(['Env-Switch: OG-OGNimbus'] ,[ 'Norm=', num2str(norm(reshape(d,[],1)))])  % OG - OGNimbus
        regModelVersion = 'default'
        
        [~,~,~,Data{4},~] = session1.plotCheckerboards(Labels,Post1Early,fh,ph(1,4),AdaptLate,flip,summaryflag);%Post1 early - Adaptation_{SS}, transition 1
        d = nanmedian(Data{4}, 4);
        title(['Transition 1'] ,[ 'Norm=', num2str(norm(reshape(d,[],1)))])
        
        [~,~,~,Data{5},~] = session1.plotCheckerboards(Labels,Post2Early,fh,ph(1,5),Post1Late,flip,summaryflag); %TM post VR early - OG post late, transition 2
        d = nanmedian(Data{5}, 4);
        title(['Transition 2'] ,[ 'Norm=', num2str(norm(reshape(d,[],1)))])
        
        [~,~,~,Data{6},~] = session2Data.plotCheckerboards(Labels,OGpostPosEarly,fh,ph(1,6),PostShort,flip,summaryflag); %TM post VR early - OG post late, transition 2
        d = nanmedian(Data{6}, 4);
        title(['Transition 3: Short Split',] ,[ 'Norm=', num2str(norm(reshape(d,[],1)))])
        
 
        set(ph(:,1),'CLim',[-1 1]*1);
        set(ph(:,2:end),'YTickLabels',{},'CLim',[-1 1]*1);
        set(ph,'FontSize',8)
        pos=get(ph(1,end),'Position');
        axes(ph(1,end))
        colorbar
        set(ph(1,end),'Position',pos);
        
        set(gcf,'color','w');
        
    
end


end