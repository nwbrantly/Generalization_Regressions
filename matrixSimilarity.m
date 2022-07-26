%getting similarity matrix 
clear all;clc

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
%%

% load('A_15_AsymC4_EarlyLateAdaptation.mat')
load ATR_4_AsymC5
aux1= sign(C(:,2));
aux2=sign(C(:,3));
% aux3=find(C(:,2)<0.1 & C(:,2)>-0.1) ;
% aux1(aux3)=1;aux2(aux3)=-1;
temp=find(aux1==aux2);


% temp(aux3)=[];

ytl={'GLU','TFL','HIP','RF','VL','VM','BF', 'SEMB','SEMT','MG','LG','SOL','PER','TA'};
% ytl=newLabelPrefix;
% ytl=ytl(end:-1:1);
yt=1:14;
fs=14;

figure 
subplot(1,2,1)
imagesc((reshape(C(:,2),12,14)'))
caxis([-1 1])
set(gca,'XTick',[],'YTick',yt,'YTickLabel',ytl,'FontSize',fs)
title('Steady State')

subplot(1,2,2)
imagesc((reshape(C(:,3),12,14)'))
caxis([-1 1])
set(gca,'XTick',[],'YTick',yt,'YTickLabel',ytl,'FontSize',fs)
title('Early Post Adapt')
colormap(flipud(map))
% sam=0;
for i=1:length(aux1)
    if sum(i==temp)==1
%         idx=find(i==temp);
        similar_matrix(i,1)=C(i,2);
%         sam=sam+1;
        diff_matrix1(i,1)=0;
        diff_matrix2(i,1)=0;
    else
        similar_matrix(i,1)=0;
        diff_matrix1(i,1)=C(i,2);
        diff_matrix2(i,1)=C(i,3);
    end
end

Cnew=[similar_matrix diff_matrix1 diff_matrix2];


% 
figure
subplot(1,3,1)
imagesc((reshape(Cnew(:,1),12,14)'))
caxis([-1 1])
set(gca,'XTick',[],'YTick',yt,'YTickLabel',ytl,'FontSize',fs)
title(['Similar Matrix'],[ '%Similary=' num2str((length(temp) *100)/168)])


subplot(1,3,2)
imagesc((reshape(Cnew(:,2),12,14)'))
caxis([-1 1])
set(gca,'XTick',[],'YTick',yt,'YTickLabel',ytl,'FontSize',fs)
title('Diff 1')

subplot(1,3,3)
imagesc((reshape(Cnew(:,3),12,14)'))
caxis([-1 1])
set(gca,'XTick',[],'YTick',yt,'YTickLabel',ytl,'FontSize',fs)
title('Diff 2')


colormap(flipud(map))

%%

%% Free model - Linear regression - Asymmetry 
 

% load('ATS_11_Asym_EarlyLateAdaptation.mat')
% load('allDataRedAlt_ATS_fixCandD1_280322T212228.mat')
% load('AUFV1_5_Asym_EarlyLateAdaptation')
% load('ATS_11_Asym_EarlyLate5Adaptation.mat')
% load('AUFV4_5_Asym_EarlyLateAdaptation.mat')
% load('NTS_5_Asym_AdaptationPosShort.mat')
% load('NTR_4_Asym_AdaptationPosShort.mat')
% load('ATR_4_Asym_EarlyLate10Adaptation.mat')

% load('CTS_5_Asym_AdaptationPosShort.mat')
% load('CTR_4_Asym_AdaptationPosShort.mat')
% load PATR_2_AsymC3_EarlyLateAdaptation
% figure 
% fname='dynamicsData_PATR.h5';
%  fname= 'dynamicsData_C_s12V2.h5';
fname='dynamicsData_ATR_V4.h5';
EMGdata=h5read(fname,'/EMGdata');
 
binwith=1;
[Y,Yasym,Ycom,U,Ubreaks]=groupDataToMatrixForm(1:size(EMGdata,3),0,fname);
Uf=[U;ones(size(U))];
% C=[C1 C3];
C=Cnew;
% C=[C(:,1) C(:,2)];

model.C=C;
% Y=datSet.out;
Cinv=pinv(model.C)';
X2asym = Yasym*Cinv; %x= y/C
Y2asym= C * X2asym' ; %yhat = C 
   
figure
subplot(3,1,1)
hold on
plot( movmean(X2asym(:,1),binwith))
% pp=patch([50 600 600 50],[-0.5 -0.5 1.6 1.6],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
% uistack(pp,'bottom')


% title('C_1 = Early Adaptation')
% yline(0)

subplot(3,1,2)
plot(movmean(X2asym(:,2),binwith))

subplot(3,1,3)
plot(movmean(X2asym(:,3),binwith))
hold on
legend('reactive','adaptive','AutoUpdate','off')
% pp=patch([50 939 939 50],[-0.5 -0.5 1.6 1.6],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
% pp=patch([50 600 600 50],[-0.5 -0.5 1.6 1.6],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
% uistack(pp,'bottom')
axis tight
% yline(0)
% yline(1)
% title('C_2 = Late Adaptation')

set(gcf,'color','w')

% subplot(2,1,2)
% RMSEasym= sqrt(mean((Yasym-Y2asym').^2,2));
% plot(movmean(RMSEasym,binwith))
% ylabel('RMES')
% % pp=patch([50 939 939 50],[0 0 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
% pp=patch([50 600 600 50],[-0.5 -0.5 1.6 1.6],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
% uistack(pp,'bottom')
% axis tight
% end
% legacy_vizSingleModelMLMC_FreeModel(model,Yasym',Uf)
% load ATR_4_AsymC5
% model.C=C(:,1:2);
% legacy_vizSingleModelMLMC_FreeModel(model,Yasym',Uf)