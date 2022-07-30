% Ploting data used for the R01 

% load('allDataRedAlt_ATS_fixCandD1_280322T212228.mat')
figure
labels={'YA_{TR} n=4','YA_{TS} n=11','OA_{TR} n=13','OA_{TS} n=5',...
    'Stroke_{TR} n=14','OAS4_{TS} n=5','OAS4_CV2_{TS} n=5 '};
l=0;
group=[3];
for g=group
    l=l+1;
    xhat=[];
%     x=[];
    C=[];
    Y=[];
    if g==1
        %         load('allDataRedAlt_PATS_fixCandD_Adaptation040422T154503.mat')
        load('YA_TR_fixDandCV4_20220316T114557.mat')
        color=[0 0.4470 0.7410];
           color2=[0.8500 0.3250 0.0980];
           color3=[0.9290 0.6940 0.1250];
        
    elseif g==2
        %         load('allDataRedAlt_PATS_fixCandD_Perturbations040422T154842.mat')
        load('YA_TS_fixCandD1_280322T212228.mat')
        color=[0.8500 0.3250 0.0980];
        
    elseif g==3
        load('OA_TR_fixDandC_160322T155119.mat')
        color=[0.9290 0.6940 0.1250];
    elseif g==4
        load('OA_TS_fixCandD_Adaptation070422T142847.mat')
        color=[0.4940 0.1840 0.5560];
    elseif g==5
        load('Stroke_TR_fixCandD_040422T114934.mat')
        color= [0.4660 0.6740 0.1880];
    elseif g==6
        load('OAV4_TS_fixCandD_Adaptation070422T162112.mat')
        color=[0.3010 0.7450 0.9330];
    elseif g==7
        load('OAV4_TS_fixCV2andD_Adaptation070422T162112.mat')
        color=[0.6350 0.0780 0.1840];
        
    end
    ind=find(diff(datSet.in(1,:))~=0);
    
%     C3=nanmean(datSet.out(:,1:40),2);
%     C=[modelRed.C];
    Cinv=pinv(C)';
    Y=datSet.out;
%         Y=Y(:,1:1150);
    %     Y(:,ind)=nan;
    xhat = Y'*Cinv;
    
for step=1:length(Y)
    xhat2(step,:) = lsqr(C,Y(:,step));
%     x(step,:) = lsqlin(C,Y(:,step),C,Y(:,step)) ;
%     x(:,step) = lsqlin(C,datSet.out(:,step),C,datSet.out(:,step)) ;
end
    %     figure
 %%   
    subplot(2,1,1)
    plot( movmean(xhat(:,1),5),'LineWidth',2,'Color',color)
    hold on
%     plot( movmean(xhat2(:,1),5),'--','LineWidth',2,'Color',color2)
%     plot( movmean(x(1,:),5),'LineWidth',2,'Color',color3)
    hold on
    ylabel('X_{reactive}')
    xlabel('strides')
    title('C_1 = Early Adaptation')
    pp=patch([50 950 950 50],[-1 -1 2 2],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
    uistack(pp,'bottom') 
    yline(0)
    
    subplot(2,1,2)
    hold on
    plot(movmean(xhat(:,2),5),'LineWidth',2,'Color',color,'DisplayName',labels{g})
%         plot( movmean(xhat2(:,2),5),'--','LineWidth',2,'Color',color2)
%          plot( movmean(x(2,:),5),'LineWidth',2,'Color',color3)
%     legend('Location','NorthEastOutside')
    
    ylabel('X_{learning}')
    xlabel('strides')
    yline(0.1)
    title('C_2 = Late Adaptation')
    pp=patch([50 950 950 50],[-.5 -0.5 2 2],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none','DisplayName','Adaptation');
    uistack(pp,'bottom') 
    
%     
%      subplot(3,1,3)
%      hold on
%     li{l}=plot(movmean(xhat(:,3),5),'LineWidth',2,'Color',color,'DisplayName',labels{g})
%         plot( movmean(xhat2(:,3),5),'--','LineWidth',2,'Color',color2)
%     legend('Location','NorthEastOutside')
    
%     ylabel('X_{tied}')
%     xlabel('strides')
%     yline(1,'DisplayName','Baseline')
%     title('C_3 =Baseline ')
%     pp=patch([50 950 950 50],[0 0 2 2],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none','DisplayName','Adaptation');
    uistack(pp,'bottom') 
    
    
end

% legend([li{:}]',[labels])
% legend('OA_{TR}','OA_{TS}')
% legend('YA_{TS}','OA_{TS}')

set(gcf,'color','w')
%%
% load('allDataRedAlt_ATS_fixCandD1_280322T212228.mat')
% ind=find(diff(datSet.in(1,:))~=0);
% Y(:,ind)=nan;
C=modelRed.C;
Cinv=pinv(C)';
Y=datSet.out;
% Y(:,ind)=nan;
xhat = Y'*Cinv;
xhat_ats=xhat;

legacy_vizSingleModelMLMC_FreeModel(modelRed,datSet.out,datSet.in)

load('allDataRedAlt_fixDandCV4_20220316T114557.mat')
%
% ind=find(diff(datSet.in(1,:))~=0);
% Y(:,ind)=nan;
C=modelRed.C;
Cinv=pinv(C)';
Y=datSet.out;
% Y(:,ind)=nan;
xhat = Y'*Cinv;
legacy_vizSingleModelMLMC_FreeModel(modelRed,datSet.out,datSet.in)

%%
% figure
subplot(2,3,1)
plot( movmean(xhat(:,1),5))
hold on
plot(xhat_ats(:,1))
yline(0)
title('C_1 = Early Adaptation')
legend('ATR','ATS')

subplot(2,3,2)
plot(movmean(xhat(:,2),1))
hold on
plot(movmean(xhat_ats(:,2),5))
yline(0)
title('C_2 = Late Adaptation')
legend('ATR','ATS')

set(gcf,'color','w')

%% R - Squared 

model{1}.C=modelRed.C;
Y2=datSet.out;
C=model{1}.C;
Cinv=pinv(C)';
xhat = Y2'*Cinv; %x= y/C
yhat= C * xhat' ; %yhat = C 
for step= 1:size(Y2,2)
%     RMSE(step,1) = sqrt(sum((Y2(:,step) - yhat(:,step)).^2)/size(Y2,1))
%     RMSE2(step,1) = sqrt(sum((Y2(:,step) - yhat(:,step)).^2))
%     RMSE3(step,1) = sqrt(mean(((Y2(:,step) - yhat(:,step)).^2)))
    RMSE(step,1) = sqrt(immse(Y2(:,step),yhat(:,step)));
    Rsq(step,1) = 1 - sum((Y2(:,step) - yhat(:,step)).^2)/sum((Y2(:,step) - mean(Y2(:,step))).^2);
end



  Rsq = 1 - sum((Y2- yhat).^2)/sum((Y2- mean(Y2)).^2);

%%
figure
speed=[zeros(50,1) ; ones(30,1); zeros(50,1); -1*ones(30,1) ; zeros(100,1); linspace(0,1,10)'; ones(30,1); zeros(50,1) ;...
    ones(450,1);zeros(300,1); ones(30,1);zeros(20,1);...
    rand(1)*ones(30,1);zeros(20,1);rand(1)*ones(30,1);zeros(20,1);rand(1)*ones(30,1);zeros(20,1);rand(1)*ones(30,1);zeros(20,1);...
    ones(30,1);zeros(20,1);rand(1)*ones(30,1);zeros(20,1);rand(1)*ones(30,1);zeros(20,1);ones(30,1);zeros(20,1);...
    rand(1)*ones(30,1);zeros(20,1);ones(150,1);zeros(100,1);];%...
%     -1*[ones(450,1);zeros(300,1); ones(30,1);zeros(20,1);...
%     ones(30,1);zeros(20,1);ones(30,1);zeros(20,1);ones(30,1);zeros(20,1);ones(30,1);zeros(20,1);...
%     ones(30,1);zeros(20,1);ones(30,1);zeros(20,1);ones(30,1);zeros(20,1);ones(30,1);zeros(20,1);...
%     ones(30,1);zeros(20,1);ones(150,1);zeros(100,1)]];

plot(speed,'LineWidth',2,'Color','k')
ylabel({'speed difference';'m/s'})
xlabel('strides')
set(gcf,'color','w')
ylim([-.8 .9])
% axis([0 1800 -1.1 1.1])

%%


