%%
% upluad your path 
addpath(genpath('/Users/dulcemariscal/Documents/GitHub/Generalization_Regressions'))
addpath(genpath('/Users/dulcemariscal/Documents/GitHub/labTools'))
addpath(genpath('/Users/dulcemariscal/Documents/GitHub/LongAdaptation'))
addpath(genpath('/Users/dulcemariscal/Documents/GitHub/R01'))
addpath(genpath('/Users/dulcemzariscal/Documents/GitHub/splitbelt-EMG-adaptation'))
addpath(genpath('/Users/dulcemariscal/Documents/GitHub/EMG-LTI-SSM'))
addpath(genpath('/Users/dulcemariscal/Documents/GitHub/matlab-linsys'))
rmpath(genpath('/Users/dulcezmariscal/Documents/GitHub/PittSMLlab'))
 

%%

%% Load real data:
% sqrtFlag=false;
% % subjIdx=[2:6,8,10:15]; %Excluding C01 (outlier), C07, C09 (less than 600 strides of Post), C16 (missed first trial of Adapt)
% % for s=[1:4]
%     
% subjIdx=[1:5];   
% % subjIdx=[];%
% [Y,Yasym,Ycom,U,Ubreaks]=groupDataToMatrixForm(subjIdx,sqrtFlag);
% Uf=[U;ones(size(U))];
% % Looking in to asymmetry 
% datSet=dset(Uf,Yasym');

%looking into individual legs
% datSet=dset(Uf,Y');
% binwith=5;
%% Reduce data
% Y=datSet.out;
% U=datSet.in;
% X=Y-(Y/U)*U; %Projection over input
% s=var(X'); %Estimate of variance
% flatIdx=s<.005; %Variables are split roughly in half at this threshold

clear all;clc;close all
%% Free model - Linear regression - Asymmetry with baseline

% %Testing 
% fname='dynamicsData_BATR_subj_12_RemoveBadMuscles1_splits_0_V4.h5';
% fname='dynamicsData_BAT_subj_24_RemoveBadMuscles1_splits_0_V4.h5'
fname='dynamicsData_BATR_subj_12_RemoveBadMuscles1_splits_0_WithPost2V2.h5'
% fname='dynamicsData_BATS_subj_12_RemoveBadMuscles1_splits_0_WithPost2V2_WogBaseline.h5'

%posterior muscles
% fname='dynamicsData_BATR_subj_12_RemoveBadMuscles1_splits_0_PosteriorMuscles.h5';
% load BAT_24_IndvLegsC17_ShortPertubations_RemovedBadMuscle_1RemoveBias_0_PosteriorMuscles.mat

%BATR 
% load BATR_12_AsymC16_ShortPertubations_RemovedBadMuscle_1.mat
% load BATR_12_IndvLegsC16_ShortPertubations_RemovedBadMuscle_1RemovBias_0.mat

%BATS
% load BATS_12_IndvLegsC17_ShortPertubations_RemovedBadMuscle_1RemovBias_0.mat

% ALL 24 participatns
% load BAT_24_AsymC16_ShortPertubations_RemovedBadMuscle_1.mat
load BAT_24_IndvLegsC16_ShortPertubations_RemovedBadMuscle_1.mat


EMGdata=h5read(fname,'/EMGdata');
 
binwith=5;
[Y,Yasym,~,U,~,Ysum,Yinv]=groupDataToMatrixForm_Update(1:size(EMGdata,3),fname,0);
Yasym=Yinv;
% Yasym=Yasym(:,1:size(Yasym,2)/2); %SLOW
% Yasym=Yasym(:,169:end); %FAST

range=1:12;
% Yasym=Yasym(:,range); %One mucles
% Yasym=Yinv(1:480,:);
% Yasym=Y;
% [Y,Yasym,Ycom,U,Ubreaks,Ysum]
% Yasym=Yasym(481:end,:);
% U=U(481:end);
Uf=[U;ones(size(U))];
% Uf=Uf(:,1:480);
% Uf=Uf(:,41:end);
% Y=datSet.out;
%  
% U=datSet.in;
% Code to remove muscle with low variability 
% Uinv=pinv(U)';
% X=Yasym'-(Yasym'*Uinv')*U; %Projection over input
% s=var(X'); %Estimate of variance
% flatIdx=s<.005; %Variables are split roughly in half at this threshold
% Yasym=Yasym(:,~flatIdx);


% C=[C1 C3];
% C=Cnew;
reactive=find(strcmp(epochOfInterest,'Ramp')==1);
reactive2=find(strcmp(epochOfInterest,'Tied post ramp')==1);
context= find(strcmp(epochOfInterest,'Optimal')==1);
base= find(strcmp(epochOfInterest,'TM base')==1);
% base= find(strcmp(epochOfInterest,'OG base')==1);
% Casym=[C(:,reactive) C(:,context) -C(:,reactive2) C(:,base)]; % EMGreactive and EMGcontext

% reactive2=find(strcmp(epochOfInterest,'NegShort_{late}')==1);
% reactive=find(strcmp(epochOfInterest,'Adaptation_{early}')==1);
% context= find(strcmp(epochOfInterest,'Adaptation')==1);

% reactive=find(strcmp(epochOfInterest,'Post1_{Early}')==1);
% context= find(strcmp(epochOfInterest,'Adaptation')==1);


NNMF=0
PCA_analysis=0
removemean=0
const=0
removebaseline=1


% mix=[C(1:168,context); C(169:end,reactive2)];
Casym=[C(:,reactive2) C(:,context)]; % EMGreactive and EMGcontext
% Casym=[C(:,context)];
% Cnorm=vecnorm([C(:,reactive2) C(:,context)]);
% % Casym=[C(:,reactive) C(:,context) -C(:,reactive2)]; % EMGreactive and EMGcontext
% Casym= [C(:,base)]; 
% Casym=[ C(:,context) -C(:,reactive2)]; % EMGreactive and EMGcontext

% Casym=[C(1:size(C,1)/2,reactive2)  C(1:size(C,1)/2,context)]; %SLOW
% Casym=[C(169:end,reactive2)  C(169:end,context)]; %FAST

% Casym=[C(range,reactive2)  C(range,context)]; %one muscle
% Casym=[ C(range,context)]; %one muscle
% 
% load BATR_12_IndvLegsC13_ShortPertubations_RemovedBadMuscle_0.mat
% Csum=C(:,1)+fftshift(C(:,1),2);
% Csum=Csum(1:size(Csum,1)/2,1);
% 
% %getting the norm 
% EMGsumNorm=norm(Csum)
% for s=1:length(Yasym)
% data(s)=norm(Yasym(s,:));
% end
% EMGsumNorm=norm(Csum)./mean(data);
% Yasym=Yasym(41:end,:);

if const==1
    Casym=[Casym Csum./mean(data)];
    Ymodel=[Csum'./mean(data)+ Yasym]' ;
else
%     C=[Casym];
    Ymodel=Yasym';
    
end

if removebaseline==1
    bias=nanmean(Yasym(5:30,:));
    Casym=Casym-bias';
%     Casym=Casym./vecnorm(Casym);
    Ymodel=Ymodel-bias'; % Transpose the EMG data to match equations
%      load BAT_24_IndvLegsC16_ShortPertubations_RemovedBadMuscle_1.mat
%      model.C1=[C(:,reactive) C(:,context)]-bias';
%      model.C2=[C(:,reactive2) C(:,context)]-bias';
    
end

if removemean==1
    m=nanmean(Yasym(41:end,:),1); %m stands fro mean
    Casym=Casym-m';
    Ymodel=Ymodel-m'; % Transpose the EMG data to match equations

end


if PCA_analysis==1
    [pp,cc,aa]=pca(Ymodel(:,35:485)');%,'Centered','off');
    PC= [(cc(:,1:2)*pp(:,1:2)') + nanmean(Ymodel')]';
    X=cc(:,1:2);
    C=pp(:,1:2);
    
    
    model.C=C;
    X2asym =X;
    
    
elseif NNMF==1
    
    swift=abs(min(Ymodel,[],'all'));
    data=Ymodel+swift;
    [W,H] = nnmf(data',2);
    data_hat=W*H;
    data_hat2=data_hat-swift;
    
    model.C=H';
    X2asym =W;
     
    
    
else
    model.C=Casym;
    Cinv=pinv(model.C);
    X2asym = Cinv*Ymodel; %x= y/C
    Y2asym=  model.C * X2asym  ; %yhat = C
    

    X2asym=[X2asym'];
end

% Xasym=[X2asym(1:40,:); nan(1,size(X2asym,2));X2asym(41:480,:);nan(1,size(X2asym,2));X2asym(481:end,:);nan(1,size(X2asym,2));X2asym(681:end,:)];
% Xasym=[X2asym(1:40,:); nan(1,size(X2asym,2));X2asym(41:480,:);nan(1,size(X2asym,2));X2asym(481:end,:)];

% Xasym=[X2asym(1:200,:);nan(1,size(X2asym,2));X2asym(201:end,:)];

% Xasym=[X2asym(1:40,:); nan(1,2);X2asym(41:480,:);nan(1,2);X2asym(480:end,:)];
fdr=0.1;
% [pvalc,hc,alphaAdj_c]=checkerstatsV2(reshape(C(:,1),12,14)',[],1,0,fdr,'benhoch',0);
%%

% Xasym=[X2asym(41:480,:);nan(1,size(X2asym,2));X2asym(481:681,:);nan(1,size(X2asym,2));X2asym(681:end,:)];
Xasym=[X2asym(1:40,:);nan(1,size(X2asym,2));X2asym(41:480,:);nan(1,size(X2asym,2));X2asym(481:end,:)];
% Xasym=[X2asym(1:40,:);nan(1,size(X2asym,2));X2asym(41:80,:);...
%     nan(1,size(X2asym,2));X2asym(81:520,:);nan(1,size(X2asym,2));X2asym(521:end,:)];
% Xasym=[X2asym(681:843,:);nan(1,size(X2asym,2))];

figure
subplot(4,1,1)
hold on
% scatter(1:length(movmean(Xasym(:,1),binwith)), movmean(Xasym(:,1),binwith),10,'k','filled')
scatter(1:length(movmean(Xasym(:,1),binwith)), movmean(Xasym(:,1),binwith),'filled','MarkerFaceColor',"#EDB120") %"#77AC30" )%
pp=patch([40 480 480 40],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
% plot( movmean(Xasym(:,1),binwith),'Color',"#EDB120",'LineWidth',5)
% pp=patch([40 480 480 40],[0.5 0.5 1.5 1.5],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
% pp=patch([40 480 480 40],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
% pp=patch([0 440 440 0],[-0.5 -0.5 1 1],.7*ones(1reactive2,3),'FaceAlpha',.2,'EdgeColor','none');

% pp=patch([40 480 480 40],[-4 -4 5 5],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');

% legend('Baseline','AutoUpdate','off')
legend('Reactive','AutoUpdate','off')
%     legend('Removal Perturbation','AutoUpdate','off')
% uistack(pp,'bottom')
yline(0)
ylabel({'Reactive';'(A.U)'})
%     ylabel({'Removal';'(A.U)'})
xlabel('strides')

if size(X2asym,2)>=2
    % figure
    subplot(4,1,2)
    hold on
    scatter(1:length(movmean(Xasym(:,2),binwith)), movmean(Xasym(:,2),binwith),'filled','MarkerFaceColor',"#77AC30")
    pp=patch([40 480 480 40],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
    % % plot( movmean(Xasym(:,2),binwith),'Color',"#77AC30",'LineWidth',5)
    % % pp=patch([40 480 480 40],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
    % % pp=patch([0 440 440 0],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
    %
    legend('Contextual','AutoUpdate','off')
    % % legend('Switch','AutoUpdate','off')
    % % uistack(pp,'bottom')
    yline(0)
    ylabel({'Contextual';'(A.U)'})
    xlabel('strides')
end

if size(X2asym,2)>=3
    subplot(4,1,3)
    scatter(1:length(movmean(Xasym(:,3),binwith)),movmean(Xasym(:,3),binwith),'filled','MarkerFaceColor',"#0072BD")
    hold on
    %     scatter(movmean(Xasym(:,3),binwith),'Color',"#0072BD",'filled')
    hold on
    %     legend('Contextual','AutoUpdate','off')
    legend('Removal Perturbation','AutoUpdate','off')
    %     pp=patch([40 480 480 40],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
    pp=patch([40 480 480 40],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
    
    uistack(pp,'bottom')
    axis tight
    yline(0)
    ylabel({'Removal';'(A.U)'})
    xlabel('strides')
end
if size(X2asym,2)==4
    subplot(4,1,4)
    scatter(1:length(movmean(Xasym(:,4),binwith)),movmean(Xasym(:,4),binwith),'filled','MarkerFaceColor','k')
    hold on
    %     plot(movmean(Xasym(:,4),binwith))
    hold on
    legend('Baseline','AutoUpdate','off')
    pp=patch([40 480 480 40],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
    uistack(pp,'bottom')
    ylabel({'Baseline';'(A.U)'})
    axis tight
    yline(0)
end
set(gcf,'color','w')
%% Define the type of analysis that you want to do 
%0 linear regression 
%1 Using hypothesize patterns as constants
%2 PCA 
%3 removing the mean for the data (varaince reconstrucitons)
analysis=0
    
 

% legacy_vizSingleModel_FreeModel_ShortAdaptation(model,Ymodel,Uf,analysis)
legacy_vizSingleModel_FreeModel_ShortAdaptation_IndvLeg(model,Ymodel,Uf,analysis,[],[])
% legacy_vizSingxleModel_FreeModel_›ShortAdae ptation(model,Ymodel',Uf,1)
%     
% model.C=C(:,1); 
% legacy_vizSingleModel_FreeeModel_ShortAdaptation(model,Ymodel',Uf,1)

% legacy_vizSingleModel_FreeModel_ShortAdaptation_IndvLeg_cond(model,Ymodel,Uf,analysis)

%% Individual muscle reconstruction - RMSE


mm= 0:12:336;
mm2=1:12:336;
for m=1:28
    muscle(m,:)=sqrt(sum((Ymodel(mm2(m):mm(m+1),:)-Y2asym(mm2(m):mm(m+1),:)).^2));
    muscles(m,:)= 1 - sum((Ymodel(mm2(m):mm(m+1),:)-Y2asym(mm2(m):mm(m+1),:)).^2)./sum((Ymodel(mm2(m):mm(m+1),:)- mean(Ymodel(mm2(m):mm(m+1),:))).^2);
    
end

% aux3= 1 - sum((residual_2).^2)./sum((Yasym2- mean(Yasym2)).^2);
% aux3=conv(aux3,ones(1,binw)/binw,'valid'); %Smoothing
% plot(aux3,'-','LineWidth',2,'DisplayName','2 States','Color','b') ;
% ylabel({'R^{2}'})


figure(1)
muscleOrder={'GLU','TFL','HIP','RF','VL','VM','BF', 'SEMB','SEMT','MG','LG','SOL','PER','TA'};
ytl= defineMuscleListV2(muscleOrder); %List of muscle 
% ytl={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT', 'VM', 'VL', 'RF', 'HIP','TFL', 'GLU'};
ytl(end:-1:1) = ytl(:);
for m=1:14
    subplot(14,1,m)
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
for m=15:28
    subplot(14,1,m-14)
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



%% PCA residuals
residual = Ymodel-Y2asym;

[pp,cc,aa]=pca(residual','Centered','off');
% [coeff,score,latent,tsquared,explained,mu]=pca(Yasym,'Centered','off');
%%Input has to be row observation  and columns variables 
%%pp - data projected in the principal component 
%%cc - Corresponding matrix of eigenvectors
%%aa - vector of eigent values 

%Recontruction of the data 
Xcentered = pp*cc'; 
% Ymean=nanmean(Yasym,2);
Ynew=Xcentered;%+Ymean;

figure 
subplot(2,1,1)
plot(aa,'LineWidth',2)
ylabel('Eigenvalues')
xline(2,'r')
xlabel('Number of Components')
% ylim([-.5 7])
 axis tight

% yline(0.85,'r')
% xline(0.85,'r')

subplot(2,1,2)
csum=cumsum(aa);
variance_explained = csum / sum(aa);
plot(variance_explained,'LineWidth',2)
xline(2,'r')
ylabel('Variance explained')
xlabel('Number of Components')
set(gcf,'color','w')


Cres=[pp(:,1)];

% Cres=[pp(:,1)+pp(:,3) pp(:,2) pp(:,4) ];
Cinv=pinv(Cres);
model.C=Cres;
Xresidual= residual'*Cinv';

figure
scatter(1:length(movmean(Xresidual(:,1),binwith)),movmean(Xresidual(:,1),binwith),10,'k','filled')
% hold on 
% scatter(1:length(movmean(Xresidual(:,2),binwith)),movmean(Xresidual(:,2),binwith),10,'b','filled')

model.C=Cres;
legacy_vizSingleModel_FreeModel_ShortAdaptation(model,residual,Uf,0)
%% Color definition 
ex1=[1,0,0];
ex2=[0,0,1];
cc=[0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];
ex1=cc(2,:);
ex2=cc(5,:);
mid=ones(1,3);
N=100;
gamma=1.5; %gamma > 1 expands the white (mid) part of the map, 'hiding' low values. Gamma<1 does the opposite
gamma=1;
map=[flipud(mid+ (ex1-mid).*([1:N]'/N).^gamma); mid; (mid+ (ex2-mid).*([1:N]'/N).^gamma)];

ytl={'GLU','TFL','HIP','RF','VL','VM','BF', 'SEMB','SEMT','MG','LG','SOL','PER','TA'};
% ytl=newLabelPrefix;
% ytl=ytl(end:-1:1);
yt=1:14;
fs=14;

% figure 
norm(pp(:,1))
subplot(1,2,2)
imagesc((reshape(pp(:,1),12,14)'))
caxis([-.4 .4])
set(gca,'XTick',[],'YTick',yt,'YTickLabel',ytl,'FontSize',fs)
title('PC1- 2 states')
colormap(map)

%%
%% Free model - Linear regression - Asymmetry 
% load PATS_3_AsymC4_ShortPertubations
% load PATR_4_AsymC4_ShortPertubations
% load PATS_3_AsymC5_ShortPertubations_RemovedBadMuscle_1
% load ATR_4_AsymC5
% figure 

% fname='dynamicsData_PATR_subjects_4.h5';
% fname='dynamicsData_PATS_subj_3_RemoveBadMuscles1.h5'
% fname='dynamicsData_PATS_subjects_3.h5';
% fname= 'dynamicsData_PATR_subj_4_RemoveBadMuscles0_splits_1.h5';
% fname='dynamicsData_ATR_V4.h5'


fname='dynamicsData_BATR_subj_4_RemoveBadMuscles0_splits_0.h5';
load BATR_4_IndvLegsC8_ShortPertubations_RemovedBadMuscle_0.mat 

EMGdata=h5read(fname,'/EMGdata');
 
binwith=5;
[Y,Yasym,~,U]=groupDataToMatrixForm(1:size(EMGdata,3),0,fname);
Uf=[U;ones(size(U))];
% C=[C1 C3];
% C=Cnew;
bias=nanmean(Yasym(5:30,:));
C=[C(:,3) C(:,4)];
C=C-bias';
Yasym=Yasym-bias;
model.C=C;
% Y=datSet.out;
Cinv=pinv(model.C)';
X2asym = Yasym*Cinv; %x= y/C
% Y2asym= C * X2asym' ; %yhat = C 

Xasym=[X2asym(1:40,:); nan(1,2);X2asym(41:490,:);nan(1,2);X2asym(490:end,:)];

figure
subplot(2,1,1)
hold on
scatter(1:length(movmean(Xasym(:,1),binwith)), movmean(Xasym(:,1),binwith),10,'k','filled')
plot( movmean(Xasym(:,1),binwith))
pp=patch([40 490 490 40],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
legend('Reactive','AutoUpdate','off')
uistack(pp,'bottom')
yline(0)

% figure
subplot(2,1,2)
hold on
scatter(1:length(movmean(Xasym(:,2),binwith)), movmean(Xasym(:,2),binwith),10,'k','filled')
plot( movmean(Xasym(:,2),binwith))
pp=patch([40 490 490 40],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
legend('Contextual','AutoUpdate','off')
uistack(pp,'bottom')
yline(0)

set(gcf,'color','w')

% legacy_vizSingleModel_FreeModel_ShortAdaptation(model,Yasym',Uf)
%%
% PCA 
[pp,cc,aa]=pca(Yasym,'Centered','off');
C=[pp(:,1:2)];
model.C=C;
Cinv=pinv(C)';
X= Yasym*Cinv;
X=[X(1:40,:); nan(1,2);X(41:490,:);nan(1,2);X(490:end,:)];

figure
subplot(2,1,1)
hold on
scatter(1:length(movmean(X(:,1),binwith)), movmean(X(:,1),binwith),10,'k','filled')
plot(movmean(X(:,1),5),'LineWidth',2)
yline(0)
title('C_1 = PCA1')


subplot(2,1,2)
hold on
scatter(1:length(movmean(X(:,1),binwith)), movmean(X(:,2),binwith),10,'k','filled')
plot(movmean(X(:,2),5),'LineWidth',2)
yline(0)
title('C_2 = PCA2')

legacy_vizSingleModel_FreeModel_ShortAdaptation(model,Yasym',Uf)
% legacy_vizSingleModelMLMC_FreeModel(model,Yasym',Uf)
%% Free model - Linear regression - Indv Legs

fname='dynamicsData_BATS_subj_12_RemoveBadMuscles0_splits_0.h5';
load BATS_12_IndvLegsC13_ShortPertubations_RemovedBadMuscle_0.mat;


EMGdata=h5read(fname,'/EMGdata');
binwith=5;
[Y,Yasym,Ycom,U,Ubreaks,Ysum]=groupDataToMatrixForm(1:size(EMGdata,3),0,fname);

Uf=[U;ones(size(U))];

Yf2= Ysum./2  + Yasym./2 ;
Ys2= Ysum./2 -  Yasym./2 ;

C=[C(:,5) C(:,6)]; 


C=[C(1:size(C,1)/2,:) C(size(C,1)/2+1:end,:)]; % ['s_{reactive}','s_{adaptive}','f_{reactive}','f_{adaptive}']

% Slow leg regression 
Cs=[C(:,1:4)];
Cinv=pinv(Cs)';
Ys=Y(:,1:size(Y,2)/2);
Yf=Y(:,size(Y,2)/2+1:end);

%Slow side 
% Cinv=pinv(C)';
Xs = Ys*Cinv; %x= y/C
Ys_hat= Cs * Xs' ; %yhat = C   

figure 
subplot(3,1,1)
plot( movmean(Xs(:,1),binwith))
hold on
plot(movmean(Xs(:,2),binwith))
plot(movmean(Xs(:,3),binwith))
plot(movmean(Xs(:,4),binwith)) 
% % plot(movmean(Xs(:,5),binwith))
% plot(movmean(Xs(:,6),binwith))
hold on
% title('C_2 = Late Adaptation')
title('S -Individual leg analysis')
set(gcf,'color','w')
yline(0) 
% legend('s_{TMbase}','s_{reactive}','s_{adaptive}')
% legend('s_{TMbase}','s_{reactive}','s_{adaptive}','f_{reactive}')
%  legend('s_{TMbase}','s_{reactive}','s_{adaptive}','f_{reactive}','f_{adaptive}','f_{TMbase}')
% legend('s_{reactive}','s_{TMbase}','f_{reactive}','f_{TMbase}')
legend('s_{reactive}','s_{adaptive}','f_{reactive}','f_{adaptive}')
% legend('s_{reactive}','s_{adaptive}','f_{reactive}')
% legend('s_{reactive}','s_{adaptive}','s_{earlypost}','f_{reactive}','f_{adaptive}','f_{earlypost}')
%  legend('s_{reactive}','s_{earlypost}','f_{reactive}','f_{earlypost}')


 
subplot(3,1,2)
R2=1 - sum((Ys'-Ys_hat).^2)./sum((Ys'- mean(Ys2')).^2);
plot(movmean(R2,binwith))
ylabel('R^{2}')

subplot(3,1,3)
RMSE= sqrt(mean((Ys-Ys_hat').^2,2));
plot(movmean(RMSE,binwith))
ylabel('RMES')

model.C=Cs;
% legacy_vizSingleModelMLMC_FreeModel(model,Ys',Uf)
legacy_vizSingleModel_FreeModel_ShortAdaptation(model,Ys',Uf,0)

%% % Fast side regression
Cf=[ C(:,3:4) C(:,1:2)];
% Cf=[C(:,1:4)];
model.C=Cf;
Cinv=pinv(Cf)';
Xf = Yf*Cinv; %x= y/C
Yf_hat= Cf * Xf' ; %yhat = C   
figure 
subplot(3,1,1)
plot( movmean(Xf(:,1),binwith))
hold on
% title('C_1 = Early Adaptation')

% subplot(2,1,2)
plot(movmean(Xf(:,2),binwith))
plot(movmean(Xf(:,3),binwith))
plot(movmean(Xf(:,4),binwith))
% plot(movmean(Xf(:,5),binwith))
% plot(movmean(Xf(:,6),binwith))
hold on
% title('C_2 = Late Adaptation')
title('F -Individual leg analysis')
set(gcf,'color','w')
yline(0)
%  legend('s_{reactive}','s_{adaptive}','s_{TMbase}','f_{reactive}','f_{adaptive}','f_{TMbase}')
% legend('s_{reactive}','s_{TMbase}','f_{reactive}','f_{TMbase}')
legend('f_{reactive}','f_{adaptive}','s_{reactive}','s_{adaptive}')
% legend('s_{reactive}','s_{adaptive}','s_{earlypost}','f_{reactive}','f_{adaptive}','f_{earlypost}')
% legend('s_{reactive}','s_{adaptive}','s_{earlypost}','f_{reactive}','f_{adaptive}','f_{earlypost}')
%  legend('s_{reactive}','s_{earlypost}','f_{reactive}','f_{earlypost}')
% legend('f_{TMbase}','f_{reactive}','f_{adaptive}')
% legend('s_{reactive}','f_{reactive}','f_{adaptive}')
subplot(3,1,2)
R2=1 - sum((Yf'-Yf_hat).^2)./sum((Yf'- mean(Yf2')).^2);
plot(movmean(R2,binwith))
ylabel('R^{2}')

subplot(3,1,3)
RMSE=[];
hold on 
RMSE= sqrt(mean((Yf-Yf_hat').^2,2));
plot(movmean(RMSE,binwith))
ylabel('RMES')

% model.C=C;
% legacy_vizSingleModelMLMC_FreeModel(model,Yf',Uf)
legacy_vizSingleModel_FreeModel_ShortAdaptation(model,Yf',Uf,0)

%% Combining left and right regressors

% C=[C1 C2 C4]; 
% Cs=[C(1:size(C,1)/2,:) C(size(C,1)/2+1:end,1)]; % ['s_{reactive}','s_{adaptive}','s_{earlypost}']
% Cinv=pinv(Cs)';
% Ys=Y(:,1:size(Y,2)/2);
% 
% 
% %Slow side 
% % Cinv=pinv(C)';
% Xs = Ys*Cinv; %x= y/C
% Ys_hat= Cs * Xs' ; %yhat = C   
% figure 
% subplot(2,1,1)
% plot( movmean(Xs(:,1),binwith))
% hold on
% % title('C_1 = Early Adaptation')
% 
% 
% 
% % subplot(2,1,2)
% plot(movmean(Xs(:,2),binwith))
% plot(movmean(Xs(:,3),binwith))
% plot(movmean(Xs(:,4),binwith))
% hold on
% % title('C_2 = Late Adaptation')
% title('S -Individual leg analysis')
% set(gcf,'color','w')
% yline(0)
% legend('s_{reactive}','s_{adaptive}','s_{EarlyPost}')
% legend('s_{reactive}','s_{adaptive}','s_{TMbase}','f_{reactive}')
% 
% 
% subplot(2,1,2)
% RMSE= sqrt(mean((Ys-Ys_hat').^2,2));
% 
% plot(movmean(RMSE,binwith))
% ylabel('RMES')
% 
% model.C=Cs;
% legacy_vizSingleModelMLMC_FreeModel(model,Ys',Uf)
% 
% 
% Cf=[C(size(C,1)/2+1:end,:) C(1:size(C,1)/2,1)];% ['f_{reactive}','f_{adaptive}','f_{earlypost}']
% Yf=Y(:,size(Y,2)/2+1:end);
% Cinv=pinv(Cf)';
% Xf = Yf*Cinv; %x= y/C
% Yf_hat= Cf * Xf' ; %yhat = C   
% figure 
% subplot(2,1,1)
% plot( movmean(Xf(:,1),binwith))
% hold on
% % title('C_1 = Early Adaptation')
% 
% % subplot(2,1,2)
% plot(movmean(Xf(:,2),binwith))
% plot(movmean(Xf(:,3),binwith))
% plot(movmean(Xf(:,4),binwith))
% hold on
% % title('C_2 = Late Adaptation')
% title('F -Individual leg analysis')
% set(gcf,'color','w')
% yline(0)
% legend('f_{reactive}','f_{adaptive}','f_{EarlyPost}')
% legend('f_{reactive}','f_{adaptive}','f_{TMbase}','s_{reactive}')
% 
% subplot(2,1,2)
% RMSE=[];
% RMSE= sqrt(mean((Yf-Yf_hat').^2,2));
% plot(movmean(RMSE,binwith))
% 
% ylabel('RMES')
% model.C=Cf;
% legacy_vizSingleModelMLMC_FreeModel(model,Yf',Uf)
% 
