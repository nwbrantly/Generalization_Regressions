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

load PATR_4_AsymC5_ShortPertubations
% load PATS_3_AsymC5_ShortPertubations_RemovedBadMuscle_1

% figure 

fname='dynamicsData_PATR_subjects_4.h5';
% fname='dynamicsData_PATS_subjects_3.h5';
% fname='dynamicsData_PATS_subj_3_RemoveBadMuscles1.h5'
EMGdata=h5read(fname,'/EMGdata');
 
binwith=5;
[Y,Yasym,~,U]=groupDataToMatrixForm(1:size(EMGdata,3),0,fname);
Uf=[U;ones(size(U))];

% Y=datSet.out;
% U=datSet.in;
Uinv=pinv(U)';
X=Yasym'-(Yasym'*Uinv')*U; %Projection over input
s=var(X'); %Estimate of variance
flatIdx=s<.005; %Variables are split roughly in half at this threshold
% Yasym=Yasym(:,~flatIdx);
% C=[C1 C3];
% C=Cnew;
bias=0;%nanmean(Yasym(5:30,:));
C=[C(:,1) C(:,4) C(:,5)];
C=C-bias';
% C=C(~flatIdx,:);
Yasym=Yasym-bias;
model.C=C;
% Y=datSet.out;
Cinv=pinv(model.C)';
X2asym = Yasym*Cinv; %x= y/C
Y2asym= C * X2asym' ; %yhat = C 

% Xasym=[X2asym(1:40,:); nan(1,3);X2asym(41:490,:);nan(1,3);X2asym(490:end,:)];
Xasym=[X2asym(1:40,:); nan(1,3);X2asym(41:490,:);nan(1,3);X2asym(490:end,:)];

figure
subplot(3,1,1)
hold on
scatter(1:length(movmean(Xasym(:,1),binwith)), movmean(Xasym(:,1),binwith),10,'k','filled')
plot( movmean(Xasym(:,1),binwith))
pp=patch([40 490 490 40],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
legend('Baseline','AutoUpdate','off')
uistack(pp,'bottom')
yline(0)

% figure
subplot(3,1,2)
hold on
scatter(1:length(movmean(Xasym(:,2),binwith)), movmean(Xasym(:,2),binwith),10,'k','filled')
plot( movmean(Xasym(:,2),binwith))
pp=patch([40 490 490 40],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
legend('Reactive','AutoUpdate','off')
uistack(pp,'bottom')
yline(0)

subplot(3,1,3)
scatter(1:length(movmean(Xasym(:,3),binwith)),movmean(Xasym(:,3),binwith),10,'k','filled')
hold on
plot(movmean(Xasym(:,3),binwith))
hold on
legend('Contextual','AutoUpdate','off')
pp=patch([40 490 490 40],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
uistack(pp,'bottom')
axis tight
yline(0)

set(gcf,'color','w')

legacy_vizSingleModel_FreeModel_ShortAdaptation(model,Yasym',Uf)
%% Free model - Linear regression - Asymmetry 
% load PATS_3_AsymC4_ShortPertubations
load PATR_4_AsymC4_ShortPertubations
% load PATS_3_AsymC5_ShortPertubations_RemovedBadMuscle_1
% load ATR_4_AsymC5
% figure 

% fname='dynamicsData_PATR_subjects_4.h5';
% fname='dynamicsData_PATS_subj_3_RemoveBadMuscles1.h5'
% fname='dynamicsData_PATS_subjects_3.h5';
fname= 'dynamicsData_PATR_subj_4_RemoveBadMuscles0_splits_1.h5';
% fname='dynamicsData_ATR_V4.h5'
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
fname='dynamicsData_PATS_subjects_3.h5';
EMGdata=h5read(fname,'/EMGdata');
 
binwith=5;
[Y,Yasym,Ycom,U,Ubreaks]=groupDataToMatrixForm(1:size(EMGdata,3),0,fname);
Uf=[U;ones(size(U))];
% datSet=dset(Uf,Yasym');
% if ss==1 
C=[C2 C3]; 
% elseif ea==1
% 
% 
% end

C=[C(1:size(C,1)/2,:) C(size(C,1)/2+1:end,:)]; % ['s_{reactive}','s_{adaptive}','f_{reactive}','f_{adaptive}']
Cinv=pinv(C)';
Ys=Y(:,1:size(Y,2)/2);
Yf=Y(:,size(Y,2)/2+1:end);

%Slow side 
% Cinv=pinv(C)';
Xs = Ys*Cinv; %x= y/C
Ys_hat= C * Xs' ; %yhat = C   
figure 
subplot(2,1,1)
plot( movmean(Xs(:,1),binwith))
hold on
% title('C_1 = Early Adaptation')



% subplot(2,1,2)
plot(movmean(Xs(:,2),binwith))
plot(movmean(Xs(:,3),binwith))
plot(movmean(Xs(:,4),binwith))
% plot(movmean(Xs(:,5),binwith))
% plot(movmean(Xs(:,6),binwith))
hold on
% title('C_2 = Late Adaptation')
title('S -Individual leg analysis')
set(gcf,'color','w')
yline(0)
%  legend('s_{reactive}','s_{adaptive}','s_{TMbase}','f_{reactive}','f_{adaptive}','f_{TMbase}')
% legend('s_{reactive}','s_{TMbase}','f_{reactive}','f_{TMbase}')
% legend('s_{reactive}','s_{adaptive}','f_{reactive}','f_{adaptive}')
% legend('s_{reactive}','s_{adaptive}','s_{earlypost}','f_{reactive}','f_{adaptive}','f_{earlypost}')
 legend('s_{reactive}','s_{earlypost}','f_{reactive}','f_{earlypost}')

 
 
subplot(2,1,2)
RMSE= sqrt(mean((Ys-Ys_hat').^2,2));

plot(movmean(RMSE,binwith))
ylabel('RMES')

model.C=C;
legacy_vizSingleModelMLMC_FreeModel(model,Ys',Uf)

% Fast side 
Xf = Yf*Cinv; %x= y/C
Yf_hat= C * Xf' ; %yhat = C   
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
 legend('s_{reactive}','s_{earlypost}','f_{reactive}','f_{earlypost}')
% 
subplot(2,1,2)
RMSE=[];
RMSE= sqrt(mean((Yf-Yf_hat').^2,2));
plot(movmean(RMSE,binwith))

ylabel('RMES')
model.C=C;
legacy_vizSingleModelMLMC_FreeModel(model,Yf',Uf)

%% 

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

%%
% figure
% subplot(2,1,2)
% plot(movmean(RMSE,binwith))
% hold on 
% plot(movmean(RMSEasym,binwith))
% ylabel('RMES')
% legend('Indv Legs','Asym')

%% Fit Models
% maxOrder=3; %Fitting up to 10 states
% %Opts for indentification:
% opts.robustFlag=false;
% opts.outlierReject=false;
% opts.fastFlag=200; %Cannot do fast for NaN filled data, disable here to avoid a bunch of warnings.
% opts.logFlag=true;
% opts.indD=[];
% opts.indB=1;
% opts.Nreps=20;
% opts.stableA=true;
% % opts.fixR=median(s(~flatIdx))*eye(size(datSetRed.out,1)); %R proportional to eye
% % opts.fixR=corrMat(~flatIdx,~flatIdx); %Full R
% % opts.fixR=[]; %Free R
% load('PATS_FixC.mat')
% opts.fixC=C;
% % load('D1_fastBase-slowBase.mat')
% % opts.fixD=[D1 zeros(12*14,1)];
% opts.fixD=[zeros(12*13,2)];
% opts.includeOutputIdx=find(~flatIdx); 
% [modelRed]=linsys.id(datSet,2,opts);
% %% Save (to avoid recomputing in the future)
% nw=datestr(now,'ddmmyyTHHMMSS');
% % save(['../../res/allDataRedAlt_' nw '.mat'],'modelRed', 'datSet', 'opts');
% save(['allDataRedAlt_PATS_fixCandD_Adaptation' nw '.mat'],'modelRed', 'datSet', 'opts');
%% Add dummy model:
% opts1=opts;
% opts1.fixA=eye(180);
% opts1.fixC=eye(180);
% opts1.includeOutputIdx=1:180;
% opts1.stableA=false;
% opts1.Nreps=0;
% opts1.refineMaxIter=10; %Giant model, cannot afford too many iterations
% [modelDummy]=linsys.id(datSet,180,opts1);
% dummyFit=modelDummy.fit(datSet,[],'KF');
% save('../../res/allDataDummyModel.mat','modelDummy','opts1','datSet','dummyFit')
% %%
% load ../../res/allDataRedAlt_20190510T175706.mat
% %% Compare models
% for i=1:maxOrder
% modelRed{i}.name=num2str(i-1);
% end
% fittedLinsys.compare(modelRed(1:11))
% set(gcf,'Units','Normalized','OuterPosition',[.4 .7 .6 .3])
% % saveFig(gcf,'../../fig/','allDataModelsRedAltCompare',0)
% %% visualize structure
% modelRed2=cellfun(@(x) x.canonize,modelRed,'UniformOutput',false);
% datSet.vizFit(modelRed2(1:10))
% %%
% linsys.vizMany(modelRed(2:6))