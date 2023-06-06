%%
% upluad your path 
addpath(genpath('/Users/dulcemariscal/Documents/GitHub/Generalization_Regressions'))
addpath(genpath('/Users/dulcemariscal/Documents/GitHub/labTools'))
addpath(genpath('/Users/dulcemariscal/Documents/GitHub/LongAdaptation'))
addpath(genpath('/Users/dulcemariscal/Documents/GitHub/R01'))
addpath(genpath('/Users/dulcemariscal/Documents/GitHub/splitbelt-EMG-adaptation'))
addpath(genpath('/Users/dulcemariscal/Documents/GitHub/EMG-LTI-SSM'))
addpath(genpath('/Users/dulcemariscal/Documents/GitHub/matlab-linsys')) 
% rmpath(genpath('/Users/dulcemariscal/Documents/GitHub/PittSMLlab'))


%%
clear all;clc;close all
%% Free model - Linear regression - Asymmetry with baseline


fname='dynamicsData_BATS_subj_12_RemoveBadMuscles1_splits_0_V4.h5';

load BATS_12_AsymC2_ShortPertubations_RemovedBadMuscle_1.mat %Loeading EMG_reactive and EMG_adaptve


EMGdata=h5read(fname,'/EMGdata'); %to get the EMG data fro size 
binwith=5; %running average window 
[Y,~,~,U]=groupDataToMatrixForm_Update(1:size(EMGdata,3),fname,0); % to get the data 

Yasym=Y; %Data that we want to fit 
Uf=[U;ones(size(U))]; % Speed differences 

reactive=find(strcmp(epochOfInterest,'Neg_')==1);
context= find(strcmp(epochOfInterest,'Optimal')==1);


Casym=[C(:,reactive) C(:,context)]; % EMGreactive and EMGcontext

% load BATR_12_IndvLegsC13_ShortPertubations_RemovedBadMuscle_0.mat
% Csum=C(:,1)+fftshift(C(:,1),2);
% Csum=Csum(1:size(Csum,1)/2,1);
% 
% %getting the norm 
% EMGsumNorm=norm(Csum)
% for s=1:680
% data(s)=norm(Yasym(s,:));
% end
% EMGsumNorm=norm(Csum)./mean(data);
% % Yasym=Yasym(41:end,:);
% const=0
% removebaseline=1
% if const==1
%     Casym=[Casym Csum./mean(data)];
%     Ymodel=[Csum'./mean(data)+ Yasym]' ;
% else

% C=[Casym]; % EMGreactive and EMGcontext
Ymodel=Yasym'; %Transposing the data

%     
% end
removebaseline=1;
if removebaseline==1
    bias=nanmean(Yasym(5:30,:));
    Casym=Casym-bias';
    Ymodel=Ymodel-bias'; % Transpose the EMG data to match equations
end
% C=[Csum./mean(data) Casym];
% C=[Csum./5 Casym];

% C=C-bias';
% C=C(~flatIdx,:);

% Ymodel=Csum'./5+ Yasym ;
% 
model.C=Casym; %adding EMGreactive and EMGcontext to model structure
Cinv=pinv(model.C); %Gettign the pseudoinverse of the EMGreactive and EMGcontext
Wasym = Cinv*Ymodel; %x= y/C
hatEMGasym=  model.C * Wasym  ; %yhat = C 

Wasym=[Wasym'];


Xasym=[Wasym(1:40,:); nan(1,size(Wasym,2));Wasym(41:480,:);nan(1,size(Wasym,2));Wasym(481:end,:)];

%%
Xasym=[Wasym(41:480,:);nan(1,size(Wasym,2));Wasym(481:end,:)];
figure
subplot(4,1,1)
hold on
% scatter(1:length(movmean(Xasym(:,1),binwith)), movmean(Xasym(:,1),binwith),10,'k','filled')
scatter(1:length(movmean(Xasym(:,1),binwith)), movmean(Xasym(:,1),binwith),'filled','MarkerFaceColor',"#EDB120")
% plot( movmean(Xasym(:,1),binwith),'Color',"#EDB120",'LineWidth',5)
% pp=patch([40 480 480 40],[0.5 0.5 1.5 1.5],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
% pp=patch([40 480 480 40],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
pp=patch([0 440 440 0],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');

% pp=patch([40 480 480 40],[-4 -4 5 5],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');

% legend('Baseline','AutoUpdate','off') 
legend('Reactive','AutoUpdate','off')
uistack(pp,'bottom')
yline(0)
ylabel({'Reactive';'(A.U)'})
xlabel('strides')

% figure
subplot(4,1,2)
hold on
scatter(1:length(movmean(Xasym(:,2),binwith)), movmean(Xasym(:,2),binwith),'filled','MarkerFaceColor',"#77AC30")
% plot( movmean(Xasym(:,2),binwith),'Color',"#77AC30",'LineWidth',5)
% pp=patch([40 480 480 40],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
pp=patch([0 440 440 0],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');

legend('Contextual','AutoUpdate','off')
% legend('Switch','AutoUpdate','off')
uistack(pp,'bottom')
yline(0)
ylabel({'Contextual';'(A.U)'})
xlabel('strides')

if size(Wasym,2)>=3
    subplot(4,1,3)
    scatter(1:length(movmean(Xasym(:,3),binwith)),movmean(Xasym(:,3),binwith),'filled','MarkerFaceColor',"#0072BD")
    hold on
%     scatter(movmean(Xasym(:,3),binwith),'Color',"#0072BD",'filled')
    hold on
%     legend('Contextual','AutoUpdate','off')
    legend('Removal Perturbation','AutoUpdate','off')
%     pp=patch([40 480 480 40],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
pp=patch([0 440 440 0],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');

uistack(pp,'bottom')
    axis tight
    yline(0)
    ylabel({'Removal';'(A.U)'})
xlabel('strides')
end
if size(Wasym,2)==4
    subplot(4,1,4)
    scatter(1:length(movmean(Xasym(:,4),binwith)),movmean(Xasym(:,4),binwith),10,'k','filled')
    hold on
    plot(movmean(Xasym(:,4),binwith))
    hold on
    legend('Removal Perturbation','AutoUpdate','off')
    pp=patch([40 480 480 40],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
    uistack(pp,'bottom')
    axis tight
    yline(0)
end
set(gcf,'color','w')
%% Plotting function 
legacy_vizSingleModel_FreeModel_ShortAdaptation(model,Ymodel,Uf,3) 


% legacy_vizSingxleModel_FreeModel_›ShortAda ptation(model,Ymodel',Uf,1)
%     
% model.C=C(:,1); 
% legacy_vizSingleModel_FreeeModel_ShortAdaptation(model,Ymodel',Uf,1)

%% Individual muscle reconstruction - RMSE


mm= 0:12:168;
mm2=1:12:168;
for m=1:14
    muscle(m,:)=sqrt(sum((Ymodel(mm2(m):mm(m+1),:)-hatEMGasym(mm2(m):mm(m+1),:)).^2));
    muscles(m,:)= 1 - sum((Ymodel(mm2(m):mm(m+1),:)-hatEMGasym(mm2(m):mm(m+1),:)).^2)./sum((Ymodel(mm2(m):mm(m+1),:)- mean(Ymodel(mm2(m):mm(m+1),:))).^2);
    
end

% aux3= 1 - sum((residual_2).^2)./sum((Yasym2- mean(Yasym2)).^2);
% aux3=conv(aux3,ones(1,binw)/binw,'valid'); %Smoothing
% plot(aux3,'-','LineWidth',2,'DisplayName','2 States','Color','b') ;
% ylabel({'R^{2}'})


figure(1)
ytl={'GLU','TFL','HIP','RF','VL','VM','BF', 'SEMB','SEMT','MG','LG','SOL','PER','TA'};
for m=1:7
    subplot(7,1,m)
    hold on
    plot(movmean(muscle(m,41:end),5))
%     scatter(1:length(muscle(m,:)),movmean(muscle(m,:),5),'filled')
    legend(ytl{m},'AutoUpdate','off');
    ylabel('RMSE')
    yline(nanmean(muscle(m,10:30)))
    
%     pp=patch([40 480 480 40],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
    pp=patch([0 440  440 0],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
 
    uistack(pp,'bottom')
end
set(gcf,'color','w')

figure(2)
for m=8:14
    subplot(7,1,m-7)
    hold on
        plot(movmean(muscle(m,41:end),5))
%     scatter(1:length(muscle(m,:)),movmean(muscle(m,:),5),'filled')
    legend(ytl{m},'AutoUpdate','off');
    ylabel('RMSE')
    yline(nanmean(muscle(m,10:30)))
%     pp=patch([40 480 480 40],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
    pp=patch([0 440  440 0],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
  
uistack(pp,'bottom')
end
set(gcf,'color','w')

%% Variance explained 
figure(3)
ytl={'GLU','TFL','HIP','RF','VL','VM','BF', 'SEMB','SEMT','MG','LG','SOL','PER','TA'};
for m=1:7
    subplot(7,1,m)
    hold on
    plot(movmean(muscles(m,41:end),5))
%     scatter(1:length(muscles(m,:)),movmean(muscles(m,:),5),'filled')
    legend(ytl{m},'AutoUpdate','off');
    ylabel('R^2')
    yline(nanmean(muscles(m,10:30)))
    
%     pp=patch([40 480 480 40],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
   
    pp=patch([0 440  440 0],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
 
uistack(pp,'bottom')
end
set(gcf,'color','w')

figure(4)
for m=8:14
    subplot(7,1,m-7)
    hold on
        plot(movmean(muscles(m,41:end),5))
%     scatter(1:length(muscles(m,:)),movmean(muscles(m,:),5),'filled')
    legend(ytl{m},'AutoUpdate','off');
    ylabel('R^2')
    yline(nanmean(muscles(m,10:30)))
%     pp=patch([40 480 480 40],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
    pp=patch([0 440  440 0],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
   
uistack(pp,'bottom')
end
set(gcf,'color','w')



