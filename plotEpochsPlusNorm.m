function plotEpochsPlusNorm(EpochsOfInteres,GroupData, Labels,plotIndSubjects,plotGroup)

if plotIndSubjects
    for i = 1:length(GroupData.ID)
        adaptDataSubject = GroupData.adaptData{1, i};
        
        fh=figure('Units','Normalized','OuterPosition',[0 0 1 1],'NumberTitle', 'off', 'Name', GroupData.ID{i});
        ph=tight_subplot(1,length(EpochsOfInteres),[.03 .005],.04,.04);
        flip=true;
        
        for eps=1:length(EpochsOfInteres)
            
            [~,~,~,Data{eps},~] = adaptDataSubject.plotCheckerboards(Labels,EpochsOfInteres{eps},fh,ph(1,eps),[],flip);
            vec_norm = norm(Data{eps});
            title({[EpochsOfInteres{eps}.Condition{1}] ['Norm=', num2str(norm(reshape(Data{eps},[],1)))]});
            
                
        end
              
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

if plotGroup
    
    fh=figure('Units','Normalized','OuterPosition',[0 0 1 1]);
    ph=tight_subplot(1,length(EpochsOfInteres),[.03 .005],.04,.04);
    flip=1;
    summFlag='nanmedian';
    
    
    for eps=1:length(EpochsOfInteres)
        
        [~,~,~,Data{eps},~] = GroupData.plotCheckerboards(Labels,EpochsOfInteres{eps},fh,ph(1,eps),[],flip,summFlag);
        Data{eps} = nanmedian(Data{eps}, 4);
        title({[EpochsOfInteres{eps}.Condition{1}] ['Norm=', num2str(norm(reshape(Data{eps},[],1)))]});
        
        
        
    end
end

set(ph(:,1),'CLim',[-1 1]*1);
set(ph(:,2:end),'YTickLabels',{},'CLim',[-1 1]*1);
set(ph,'FontSize',8)
pos=get(ph(1,end),'Position');
axes(ph(1,end))
colorbar
set(ph(1,end),'Position',pos);
set(gcf,'color','w');


end 