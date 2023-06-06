%%
% addpath(genpath('../../../EMG-LTI-SSM/'))
% addpath(genpath('../../../matlab-linsys/'))
% addpath(genpath('../../../robustCov/'))
%%
% clear all;clc;close all
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


%% Free model - Linear regression - Asymmetry with baseline
% load PATS_3_AsymC5_ShortPertubations

% load PATR_4_AsymC5_ShortPertubations
% load PATS_3_AsymC5_ShortPertubations_RemovedBadMuscle_1

% figure 
% fname='dynamicsData_PATR_subj_4_RemoveBadMuscles0_splits_0.h5';
% load PATR_4_AsymC8_ShortPertubations_RemovedBadMuscle_0.mat

%Feedback data 
%adaptation 
% fname='Adaptation_BAT_subj_7_RemoveBadMuscles0_splits_0.h5';
% load BAT_7_AsymC10_ShortPertubations_RemovedBadMuscle_0.mat

% Training 
fname='dynamicsData_BATR_subj_6_RemoveBadMuscles0_splits_0.h5';
load BATR_6_AsymC13_ShortPertubations_RemovedBadMuscle_0.mat 

% fname='dynamicsData_BATR_subj_6_RemoveBadMuscles0_splits_0NO_reactiveMuscles.h5'
% load BATR_6_AsymC13_ShortPertubations_RemovedBadMuscle_0WO_TA.mat

% fname='Post1_BATR_subj__RemoveBadMuscles0_splits_0.h5';
% load BATR_4_AsymC9_ShortPertubations_RemovedBadMuscle_0.mat 

%Data without TA
% fname= 'dynamicsData_BATR_subj_5_RemoveBadMuscles0_splits_0WO_TA.h5';
% load BATR_5_AsymC13_ShortPertubations_RemovedBadMuscle_0WO_TA.mat

% %Testing 
% fname='dynamicsData_BATS_subj_5_RemoveBadMuscles0_splits_0.h5';
% load BATS_5_AsymC13_ShortPertubations_RemovedBadMuscle_0.mat 





% load 'BATS_1_AsymC8_ShortPertubations_RemovedBadMuscle_0BATS02.mat'
% fname='dynamicsData_BATS_subj_1_RemoveBadMuscles0_splits_0BATS02.h5';

% load 'BAT_3_AsymC8_ShortPertubations_RemovedBadMuscle_0BATR03.mat'
% fname='AdaptationOnly_BAT_subj_3_RemoveBadMuscles0_splits_0.h5';

EMGdata=h5read(fname,'/EMGdata');
% fname='dynamicsData_PATR_sujects_4.h5';
% fname='dynamicsData_PATS_subjects_3.h5';
% fname='dynamicsData_PATS_subj_3_RemoveBadMuscles1.h5'
EMGdata=h5read(fname,'/EMGdata');
 
binwith=5;
[Y,Yasym,~,U,~,Ysum]=groupDataToMatrixForm(1:size(EMGdata,3),0,fname);
% [Y,Yasym,Ycom,U,Ubreaks,Ysum]
Uf=[U;ones(size(U))];
% Y=datSet.out;
%  
% U=datSet.in;
% Code to remove muscle with low variability 
Uinv=pinv(U)';
X=Yasym'-(Yasym'*Uinv')*U; %Projection over input
s=var(X'); %Estimate of variance
flatIdx=s<.05; %Variables are split roughly in half at this threshold
% Yasym=Yasym(:,~flatIdx);


% C=[C1 C3];
% C=Cnew;

% C=[C(:,5) C(:,6) C(:,9)];
% Casym=[ones(size(C(:,5))) C(:,5) C(:,6) C(:,9)];
Casym=[C(:,8) C(:,7) C(:,12)];
% Casym=[(1+.447)*C(:,6) C(:,10)];
% Casym=[(1+.7293)*C(:,5) C(:,10)];
% Casym=[C(:,8) C(:,7) C(:,12)]; %Upper bound
% load BATR_4_IndvLegsC8_ShortPertubations_RemovedBadMuscle_0.mat
load BATS_5_IndvLegsC13_ShortPertubations_RemovedBadMuscle_0.mat
Csum=C(:,1)+fftshift(C(:,1),2);
Csum=Csum(1:size(Csum,1)/2,1);

%getting the norm 
EMGsumNorm=norm(Csum)
for s=1:680
data(s)=norm(Yasym(s,:));
end
EMGsumNorm=norm(Csum)/mean(data)

const=0
removebaseline=1
if const==1
    Casym=[Casym Csum./mean(data)];
    Ymodel=[Csum'./mean(data)+ Yasym]' ;
else
    C=[Casym];
    Ymodel=Yasym';
    
end

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
model.C=Casym;
Cinv=pinv(model.C);
X2asym = Cinv*Ymodel; %x= y/C
Y2asym=  model.C * X2asym  ; %yhat = C 

 X2asym=[X2asym'];


Xasym=[X2asym(1:40,:); nan(1,size(X2asym,2));X2asym(41:480,:);nan(1,size(X2asym,2));X2asym(481:end,:)];
% Xasym=[X2asym(1:40,:); nan(1,2);X2asym(41:480,:);nan(1,2);X2asym(480:end,:)];
%%
figure
subplot(4,1,1)
hold on
% scatter(1:length(movmean(Xasym(:,1),binwith)), movmean(Xasym(:,1),binwith),10,'k','filled')
scatter(1:length(movmean(Xasym(:,1),binwith)), movmean(Xasym(:,1),binwith),'filled','MarkerFaceColor',"#EDB120")
% plot( movmean(Xasym(:,1),binwith),'Color',"#EDB120",'LineWidth',5)
% pp=patch([40 480 480 40],[0.5 0.5 1.5 1.5],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
pp=patch([40 480 480 40],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
% pp=patch([40 480 480 40],[-4 -4 5 5],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');

% legend('Baseline','AutoUpdate','off') 
legend('Reactive','AutoUpdate','off')
uistack(pp,'bottom')
yline(0)
ylabel({'Contextual';'(A.U)'})
xlabel('strides')

% figure
subplot(4,1,2)
hold on
scatter(1:length(movmean(Xasym(:,2),binwith)), movmean(Xasym(:,2),binwith),'filled','MarkerFaceColor',"#77AC30")
% plot( movmean(Xasym(:,2),binwith),'Color',"#77AC30",'LineWidth',5)
pp=patch([40 480 480 40],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
% legend('Reactive','AutoUpdate','off')
legend('Switch','AutoUpdate','off')
uistack(pp,'bottom')
yline(0)
ylabel({'Contextual';'(A.U)'})
xlabel('strides')

if size(X2asym,2)>=3
    subplot(4,1,3)
    scatter(1:length(movmean(Xasym(:,3),binwith)),movmean(Xasym(:,3),binwith),'filled','MarkerFaceColor',"#0072BD")
    hold on
%     scatter(movmean(Xasym(:,3),binwith),'Color',"#0072BD",'filled')
    hold on
%     legend('Contextual','AutoUpdate','off')
    legend('Removal Perturbation','AutoUpdate','off')
    pp=patch([40 480 480 40],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
    uistack(pp,'bottom')
    axis tight
    yline(0)
    ylabel({'Removal';'(A.U)'})
xlabel('strides')
end
if size(X2asym,2)==4
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

%%
Yhat_3pca=Ypca;
Residual3_pca=sqrt(mean((Ymodel-Yhat_3pca).^2));

%%
figure

e_2=Residual2(:,40:end)-Residual2_pca(:,40:end);
e_3=Residual3(:,40:end)-Residual3_pca(:,40:end);

hold on 
plot(movmean(Residual2(:,40:end),5),'LineWidth',2,'Color', "#0072BD")


plot(movmean(Residual3(:,40:end),5),'LineWidth',2,'Color',"#D95319")
pp=patch([0 440 440 0],[0.05 0.05 .3 .3],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
legend('2 patterns','3 patterns')
uistack(pp,'bottom')

ylabel('error')
xlabel('strides')
xlim([0 650])

%% RMES by muscle 
% muscle=[];
mm= 0:12:168;
mm2=1:12:168;
for m=1:14
    muscle3(m,:)=sqrt(mean((Ymodel(mm2(m):mm(m+1),41:51)-Yhat_3(mm2(m):mm(m+1),41:51)).^2));
end
%% RMES by muscle by model 
figure
 hold on
 x=1:5:100
for m=1:14
    
    bar(x(m),nanmean(muscle(m,:)),'FaceColor', "#0072BD")
    bar(x(m)+2,nanmean(muscle3(m,:)),'FaceColor',"#D95319")
    bar(x(m)+1,nanmean(muscle2PCA(m,:)),'FaceColor', "#EDB120")
    bar(x(m)+3,nanmean(muscle3PCA(m,:)),'FaceColor',"#7E2F8E")
    
    errorbar(x(m),nanmean(muscle(m,:)),std(muscle(m,:))/sqrt(5),'k')
    errorbar(x(m)+2,nanmean(muscle3(m,:)),std(muscle3(m,:))/sqrt(5),'k')
    errorbar(x(m)+1,nanmean(muscle2PCA(m,:)),std(muscle2PCA(m,:))/sqrt(5),'k')
    errorbar(x(m)+3,nanmean(muscle3PCA(m,:)),std(muscle3PCA(m,:))/sqrt(5),'k')
end

ytl={'GLU','TFL','HIP','RF','VL','VM','BF', 'SEMB','SEMT','MG','LG','SOL','PER','TA'};
yt=1:14;
fs=14;
set(gca,'XTick',x,'XTickLabel',ytl,'FontSize',fs)
legend('2 patterns','3 patterns','2PC','3PC')
ylabel('RMES per muscle')
% colormap(map)

%% RMES error by muscle 
mm= 0:12:168;
mm2=1:12:168;
for m=1:14
    musclee(m,:)=sqrt(mean((Ymodel(mm2(m):mm(m+1),481:485)-Yhat_2(mm2(m):mm(m+1),481:485)).^2))...
        - sqrt(mean((Ymodel(mm2(m):mm(m+1),481:485)-Yhat_2pca(mm2(m):mm(m+1),481:485)).^2));
    
    muscle3e(m,:)=sqrt(mean((Ymodel(mm2(m):mm(m+1),481:485)-Yhat_3(mm2(m):mm(m+1),481:485)).^2))...
        - sqrt(mean((Ymodel(mm2(m):mm(m+1),481:485)-Yhat_3pca(mm2(m):mm(m+1),481:485)).^2));
end

figure
 hold on
 x=1:3:60
for m=1:14
    
    bar(x(m),nanmean(musclee(m,:)),'FaceColor', "#0072BD")
     bar(x(m)+1,nanmean(muscle3e(m,:)),'FaceColor',"#D95319")
   errorbar(x(m),nanmean(musclee(m,:)),std(musclee(m,:))/sqrt(5),'k')
   
    errorbar(x(m)+1,nanmean(muscle3e(m,:)),std(muscle3e(m,:))/sqrt(5),'k')
end

ytl={'GLU','TFL','HIP','RF','VL','VM','BF', 'SEMB','SEMT','MG','LG','SOL','PER','TA'};
yt=1:14;
fs=14;
set(gca,'XTick',x,'XTickLabel',ytl,'FontSize',fs)
legend('2 patterns','3 patterns')
ylabel('error per muscle')

%%
legacy_vizSingleModel_FreeModel_ShortAdaptation(model,Ymodel,Uf,0)
% 
% legacy_vizSingleModel_FreeModel_ShortAdaptation(model,Ymodel',Uf,1)
% 
% model.C=C(:,1);
% legacy_vizSingleModel_FreeModel_ShortAdaptation(model,Ymodel',Uf,1)

%% PCA residuals
residual = Ymodel;%-Y2asym(:,482:491);

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


Cres=[pp(:,1:3)];

% Cres=[pp(:,1)+pp(:,3) pp(:,2) pp(:,4) ];
Cinv=pinv(Cres);
model.C=Cres;
Xresidual= residual'*Cinv';
Ypca= model.C *Xresidual'  ;

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

% load('ATS_11_IndvLegs_EarlyLateAdaptation.mat')

% fname='dynamicsData_ATR_V4.h5';
% fname='dynamicsData_ATS_V6.h5';
% fname='dynamicsData_PATS_subjects_3.h5';

% PATR
% fname='dynamicsData_PATR_subj_4_RemoveBadMuscles0_splits_0.h5';
% load PATR_4_IndvLegsC8_ShortPertubations_RemovedBadMuscle_0.mat


 %Feedback data 
% fname= 'dynamicsData_BATR_subj_2_RemoveBadMuscles0_splits_0.h5';
% load BATR_2_IndvLegsC8_ShortPertubations_RemovedBadMuscle_0.mat
fname='dynamicsData_BATR_subj_4_RemoveBadMuscles0_splits_0.h5';
load BATR_4_IndvLegsC8_ShortPertubations_RemovedBadMuscle_0.mat 

EMGdata=h5read(fname,'/EMGdata');

binwith=5;
[Y,Yasym,Ycom,U,Ubreaks,Ysum]=groupDataToMatrixForm(1:size(EMGdata,3),0,fname);

Uf=[U;ones(size(U))];

Yf2= Ysum./2  + Yasym./2 ;
Ys2= Ysum./2 -  Yasym./2 ;

% datSet=dset(Uf,Yasym');
% if ss==1 
C=[C(:,1) C(:,5) C(:,6)]; 
% elseif ea==1
% 
% 
% end

C=[C(1:size(C,1)/2,:) C(size(C,1)/2+1:end,:)]; % ['s_{reactive}','s_{adaptive}','f_{reactive}','f_{adaptive}']
%% Slow leg regression 
Cs=[C(:,[1:3 5])];
Cinv=pinv(Cs)';
Yf=Y(:,1:size(Y,2)/2);
Ys=Y(:,size(Y,2)/2+1:end);

%Slow side 
% Cinv=pinv(C)';
Xs = Ys*Cinv; %x= y/C
Ys_hat= Cs * Xs' ; %yhat = C   
figure 
subplot(2,1,1)
plot( movmean(Xs(:,1),binwith))
hold on
% title('C_1 = Early Adaptation')



% subplot(2,1,2)
plot(movmean(Xs(:,2),binwith))
plot(movmean(Xs(:,3),binwith))
plot(movmean(Xs(:,4),binwith))
% % plot(movmean(Xs(:,5),binwith))
% plot(movmean(Xs(:,6),binwith))
hold on
% title('C_2 = Late Adaptation')
title('S -Individual leg analysis')
set(gcf,'color','w')

% legend('s_{TMbase}','s_{reactive}','s_{adaptive}')
legend('s_{TMbase}','s_{reactive}','s_{adaptive}','f_{reactive}')
%  legend('s_{TMbase}','s_{reactive}','s_{adaptive}','f_{reactive}','f_{adaptive}','f_{TMbase}')
% legend('s_{reactive}','s_{TMbase}','f_{reactive}','f_{TMbase}')
% legend('s_{reactive}','s_{adaptive}','f_{reactive}','f_{adaptive}')
% legend('s_{reactive}','s_{adaptive}','s_{earlypost}','f_{reactive}','f_{adaptive}','f_{earlypost}')
%  legend('s_{reactive}','s_{earlypost}','f_{reactive}','f_{earlypost}')

yline(0) 
 
subplot(2,1,2)
RMSE= sqrt(mean((Ys-Ys_hat').^2,2));

plot(movmean(RMSE,binwith))
ylabel('RMES')

model.C=Cs;
% legacy_vizSingleModelMLMC_FreeModel(model,Ys',Uf)
legacy_vizSingleModel_FreeModel_ShortAdaptation(model,Ys',Uf)

% 
%% Fast side regression
Cf=[C(:,[4:6 2])];
model.C=Cf;
Cinv=pinv(Cf)';
Xf = Yf*Cinv; %x= y/C
Yf_hat= Cf * Xf' ; %yhat = C   
figure 
subplot(2,1,1)
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
% legend('s_{reactive}','s_{adaptive}','f_{reactive}','f_{adaptive}')
% legend('s_{reactive}','s_{adaptive}','s_{earlypost}','f_{reactive}','f_{adaptive}','f_{earlypost}')
% legend('s_{reactive}','s_{adaptive}','s_{earlypost}','f_{reactive}','f_{adaptive}','f_{earlypost}')
%  legend('s_{reactive}','s_{earlypost}','f_{reactive}','f_{earlypost}')
% legend('f_{TMbase}','f_{reactive}','f_{adaptive}')
legend('f_{TMbase}','f_{reactive}','f_{adaptive}','s_{reactive}')
subplot(2,1,2)
RMSE=[];
RMSE= sqrt(mean((Yf-Yf_hat').^2,2));
plot(movmean(RMSE,binwith))

ylabel('RMES')
% model.C=C;
% legacy_vizSingleModelMLMC_FreeModel(model,Yf',Uf)
legacy_vizSingleModel_FreeModel_ShortAdaptation(model,Yf',Uf)
%% Combining left and right regressors

C=[C1 C2 C4]; 
Cs=[C(1:size(C,1)/2,:) C(size(C,1)/2+1:end,1)]; % ['s_{reactive}','s_{adaptive}','s_{earlypost}']
Cinv=pinv(Cs)';
Ys=Y(:,1:size(Y,2)/2);


%Slow side 
% Cinv=pinv(C)';
Xs = Ys*Cinv; %x= y/C
Ys_hat= Cs * Xs' ; %yhat = C   
figure 
subplot(2,1,1)
plot( movmean(Xs(:,1),binwith))
hold on
% title('C_1 = Early Adaptation')



% subplot(2,1,2)
plot(movmean(Xs(:,2),binwith))
plot(movmean(Xs(:,3),binwith))
plot(movmean(Xs(:,4),binwith))
hold on
% title('C_2 = Late Adaptation')
title('S -Individual leg analysis')
set(gcf,'color','w')
yline(0)
legend('s_{reactive}','s_{adaptive}','s_{EarlyPost}')
legend('s_{reactive}','s_{adaptive}','s_{TMbase}','f_{reactive}')


subplot(2,1,2)
RMSE= sqrt(mean((Ys-Ys_hat').^2,2));

plot(movmean(RMSE,binwith))
ylabel('RMES')

model.C=Cs;
legacy_vizSingleModelMLMC_FreeModel(model,Ys',Uf)


Cf=[C(size(C,1)/2+1:end,:) C(1:size(C,1)/2,1)];% ['f_{reactive}','f_{adaptive}','f_{earlypost}']
Yf=Y(:,size(Y,2)/2+1:end);
Cinv=pinv(Cf)';
Xf = Yf*Cinv; %x= y/C
Yf_hat= Cf * Xf' ; %yhat = C   
figure 
subplot(2,1,1)
plot( movmean(Xf(:,1),binwith))
hold on
% title('C_1 = Early Adaptation')

% subplot(2,1,2)
plot(movmean(Xf(:,2),binwith))
plot(movmean(Xf(:,3),binwith))
plot(movmean(Xf(:,4),binwith))
hold on
% title('C_2 = Late Adaptation')
title('F -Individual leg analysis')
set(gcf,'color','w')
yline(0)
legend('f_{reactive}','f_{adaptive}','f_{EarlyPost}')
legend('f_{reactive}','f_{adaptive}','f_{TMbase}','s_{reactive}')

subplot(2,1,2)
RMSE=[];
RMSE= sqrt(mean((Yf-Yf_hat').^2,2));
plot(movmean(RMSE,binwith))

ylabel('RMES')
model.C=Cf;
legacy_vizSingleModelMLMC_FreeModel(model,Yf',Uf)

