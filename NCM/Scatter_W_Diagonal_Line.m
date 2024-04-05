%% Figures for the NCM presentation 

% load('NCM2023_Treadmill.mat')
% TM=X2asym;
% load('NCM2023_OG.mat')
% OG=X2asym;
% steps=481:485;

% load('musclesLabels.mat')
% load('BATR_indv_muscles.mat')
% TM=X2asym;
% load('BATS_indv_muscles.mat')
% OG=X2asym;
% steps=481:485;

load('musclesLabels.mat')
load('BATR_10242023V2.mat')
TM=W;

load('BATS_10242023V2.mat')
OG=W;
steps=41:45;

figure
hold on
for i=1:28
    if i<15
    Li{1}=scatter(-nanmean(TM(1,steps,i)),-nanmean(OG(1,steps,i)),100,"filled",'MarkerFaceColor', 'b');
    text(-nanmean(TM(1,steps,i)),-nanmean(OG(1,steps,i)),{labels(i).Data(2:end-1)})
    else
      Li{2}=scatter(-nanmean(TM(1,steps,i)),-nanmean(OG(1,steps,i)),100,"filled",'MarkerFaceColor', 'r')  ;
        text(-nanmean(TM(1,steps,i)),-nanmean(OG(1,steps,i)),{labels(i).Data(2:end-1)})
    end

end
xx=-1.5:0.1:.5;

plot(xx,xx,'k')
legend([Li{:}],['Slow';'Fast'])
ylabel({'OG_{Reactive}';'A.U'})
xlabel({'TM_{Reactive}';'A.U'})
axis square


%  xlim([-0.3 1.4])
%  ylim([-0.3 1.4])
 set(gcf,'color','w')

 
 figure
hold on
for i=1:28
    if i<15
    Li{1}=scatter(nanmean(TM(2,steps,i)),nanmean(OG(2,steps,i)),100,"filled",'MarkerFaceColor', 'b');
    text(nanmean(TM(2,steps,i)),nanmean(OG(2,steps,i)),{labels(i).Data(2:end-1)})
    else
      Li{2}=scatter(nanmean(TM(2,steps,i)),nanmean(OG(2,steps,i)),100,"filled",'MarkerFaceColor', 'r')  ;
        text(nanmean(TM(2,steps,i)),nanmean(OG(2,steps,i)),{labels(i).Data(2:end-1)})
    end
    
end
legend([Li{:}],['Slow';'Fast'])
ylabel({'OG_{Contextual}';'A.U'})
xlabel({'TM_{Contextual}';'A.U'})
xx=-.5:0.1:1;
yline(0)
%     plot(xx,xx)
 xlim([-0.5 1])
 ylim([-0.5 1])
axis square
 set(gcf,'color','w')
 
 
 %%

 binwith=5;
 analysis=0;
 trace=2;
 sz = 50;
 for i=11
     %     Xasym=[temp2(1:40,:,i);nan(1,size(temp2,2));temp2(41:480,:,i);nan(1,size(temp2,2));temp2(481:end,:,i)];
     
     % load('NCM2023_OG.mat')
     % temp2=X2asym;
      temp2=TM;
     % Xasym=[temp2(1,481:end,i)];
     Xasym=[temp2(1,41:end,i)];
     %      Xasym=[temp2(1:200,:,i)];
     % Xasym=[X2asym(1:40,:);nan(1,size(X2asym,2));X2asym(41:80,:);...
     %     nan(1,size(X2asym,2));X2asym(81:520,:);nan(1,size(X2asym,2));X2asym(521:end,:)];
     % Xasym=[X2asym(681:843,:);nan(1,size(X2asym,2))];
     
     figure
     subplot(2,1,1)
     hold on
     % scatter(1:length(movmean(Xasym(:,1),binwith)), movmean(Xasym(:,1),binwith),10,'k','filled')
     scatter(1:length(movmean(Xasym,binwith)), -movmean(Xasym,binwith),sz,'filled','MarkerFaceColor',"#EDB120") %"#77AC30" )%
     
     %     pp=patch([40 480 480 40],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
     % plot( movmean(Xasym(:,1),binwith),'Color',"#EDB120",'LineWidth',5)
     % pp=patch([40 480 480 40],[0.5 0.5 1.5 1.5],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
     % pp=patch([40 480 480 40],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
     % pp=patch([0 440 440 0],[-0.5 -0.5 1 1],.7*ones(1reactive2,3),'FaceAlpha',.2,'EdgeColor','none');
     
     % pp=patch([40 480 480 40],[-4 -4 5 5],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
     
     % legend('Baseline','AutoUpdate','off')
     legend('Negative')
     %     title(labels(i).Data)
     %     legend('Removal Perturbation','AutoUpdate','off')
     % uistack(pp,'bottom')
     yline(0)
     ylabel({'Reactive';'(A.U)'})
     %     ylabel({'Removal';'(A.U)'})
     xlabel('strides')
     
     if size(temp2(:,:,i),2)>=2
         % figure
         %        load('NCM2023_OG.mat')
         %         temp2=X2asym;
         %         Xasym=[temp2( trace,481:end,i)];
         temp2=OG;
         % Xasym=[temp2(1,481:end,i)];
         Xasym=[temp2(1,41:end,i)];
         subplot(2,1,1)
         hold on
         scatter(1:length(movmean(Xasym,binwith)), -movmean(Xasym,binwith),sz,'filled','MarkerFaceColor'," #00008B")
         %         pp=patch([40 480 480 40],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
         % % plot( movmean(Xasym(:,2),binwith),'Color',"#77AC30",'LineWidth',5)
         % % pp=patch([40 480 480 40],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
         % % pp=patch([0 440 440 0],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
         %
         legend('Contextual')
         % legend('Switch','AutoUpdate','off')
         % % uistack(pp,'bottom')
         yline(0)
         ylabel({'Contextual';'(A.U)'})
         xlabel('strides')
     end
     if i<=14
         isF=0;
     else
         isF=1;
     end
     set(gcf,'color','w')
     %     legacy_vizSingleModel_FreeModel_ShortAdaptation_IndvLeg(model{i},Ymuscles(:,:,i)',Uf,analysis,{labels(i).Data(2:end-1)},isF)
     set(gcf,'color','w')
     
     %     temp5=corrcoef(model{i}.C);
     %     correlation(i)=temp5(2);
     
 end