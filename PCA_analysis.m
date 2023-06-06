% load 'BAT_3_AsymC8_ShortPertubations_RemovedBadMuscle_0BATR03.mat'
% fname='AdaptationOnly_BAT_subj_3_RemoveBadMuscles0_splits_0.h5';
% fname='dynamicsData_BAT_subj_24_RemoveBadMuscles1_splits_0_V4.h5';
% fname='dynamicsData_BATS_subj_12_RemoveBadMuscles1_splits_0_V4.h5';
% fname='dynamicsData_BATR_subj_12_RemoveBadMuscles1_splits_0_V4.h5';
fname='dynamicsData_BATS_subj_12_RemoveBadMuscles1_splits_0_WithPost2V2.h5'
load BAT_24_AsymC16_ShortPertubations_RemovedBadMuscle_1.mat
% load BATR_12_AsymC15_ShortPertubations_RemovedBadMuscle_1
EMGdata=h5read(fname,'/EMGdata');
 
binwith=5;
[Y,Yasym,~,U]=groupDataToMatrixForm(1:size(EMGdata,3),0,fname);
Uf=[U;ones(size(U))];
%  Y=datSet.out;
%  Y=Y(:,1:1150);

% bias=nanmedian(Yasym(5:30,:));
% Yasym=Yasym-bias;
% [pp,cc,aa]=pca(Y');
% [pp,cc,aa]=pca(Y');
% postadapt=Yasym(481:end,:);
% adapt=Yasym(41:480,:);
% [pp,cc,aa]=pca(postadapt);%,'Centered','off');
% [pp_2,cc_,aa_]=pca(adapt);%,'Centered','off');
% [pp_3,cc_3,aa3]=pca(Yasym);
% [coeff,score,latent,tsquared,explained,mu]=pca(Yasym,'Centered','off');
%%Input has to be row observation  and columns variables 
%%pp - Corresponding matrix of eigenvectors
%%cc - data projected in the principal component 
%%aa - vector of eigent values
% 
% figure 
% subplot(3,1,1)
% csum=cumsum(aa3);
% variance_explained_all = csum / sum(aa3);
% plot(variance_explained_all,'LineWidth',2)
% title('All data')
% xline(2,'r')
% ylabel('Variance explained')
% xlabel('Number of Components')
% set(gcf,'color','w')
% 
% 
% subplot(3,1,2)
% csum=cumsum(aa_);
% variance_explained_adapt = csum / sum(aa_);
% plot(variance_explained_adapt,'LineWidth',2)
% xline(2,'r')
% title('Adaptation')
% ylabel('Variance explained')
% xlabel('Number of Components')
% set(gcf,'color','w')
% 
% 
% 
% subplot(3,1,3)
% csum=cumsum(aa);
% variance_explained_post = csum / sum(aa);
% plot(variance_explained_post,'LineWidth',2)
% xline(2,'r')
% title('Post-adaptation')
% ylabel('Variance explained')
% xlabel('Number of Components')
% set(gcf,'color','w')
% Yasym=Yasym(481:end,:);
Ymodel=Yasym';
swift=abs(min(Yasym',[],'all'));
%LMSE Dulce's fit 

reactive=find(strcmp(epochOfInterest,'Ramp')==1);
context= find(strcmp(epochOfInterest,'Optimal')==1);
reactive2=find(strcmp(epochOfInterest,'Tied post ramp')==1);
base= find(strcmp(epochOfInterest,'TM base')==1);
% Casym=[C(:,5) C(:,6)];%./vecnorm([C(:,5) C(:,6)]); 
Casym=[C(:,reactive) C(:,context) C(:,reactive2) C(:,base) ]; 
% Casym=[C(:,context) C(:,reactive2) ]; 
% Casym=[C(:,reactive)]; 
Cinv=pinv(Casym); %Gettign the pseudoinverse of the EMGreactive and EMGcontext
Wasym = Cinv*(Ymodel); %x= y/C
hatEMGasym=  (Casym * Wasym); %yhat = C 

C_dulce=Casym;
W_dulce=Wasym;


%LMSE Adaptation fit 
Casym=[C(:,8) C(:,7)]; 
Cinv=pinv(Casym); %Gettign the pseudoinverse of the EMGreactive and EMGcontext
Wasym = Cinv*Ymodel; %x= y/C
hatEMGasym_fit=  Casym * Wasym  ; %yhat = C 


%LMSE Post fit 
Casym=[C(:,12)]; 
Cinv=pinv(Casym); %Gettign the pseudoinverse of the EMGreactive and EMGcontext
Wasym = Cinv*Ymodel; %x= y/C
hatEMGasym_Pfit=  Casym * Wasym  ; %yhat = C 

%LMSE using negative pertubation
context= find(strcmp(epochOfInterest,'Optimal')==1);
reactive=find(strcmp(epochOfInterest,'NegShort_{late}')==1);
Casym=[C(:,reactive) C(:,context)]; 
Cinv=pinv(Casym); %Gettign the pseudoinverse of the EMGreactive and EMGcontext
Wasym = Cinv*Ymodel; %x= y/C
hatEMGasym_Neg=  Casym * Wasym  ; %yhat = C 

%LMSE using negative pertubation 3 states
reactive=find(strcmp(epochOfInterest,'Ramp')==1);
context= find(strcmp(epochOfInterest,'Optimal')==1);
reactive2=find(strcmp(epochOfInterest,'NegShort_{late}')==1);
Casym=[C(:,reactive) C(:,context) C(:,reactive2) ]; 
Cinv=pinv(Casym); %Gettign the pseudoinverse of the EMGreactive and EMGcontext
Wasym = Cinv*Ymodel; %x= y/C
hatEMGasym_Neg2=  Casym * Wasym  ; %yhat = C 

%LMSE using negative pertubation
context= find(strcmp(epochOfInterest,'Optimal')==1);
reactive=find(strcmp(epochOfInterest,'Tied post ramp')==1);
Casym=[C(:,reactive) C(:,context)]; 
Cinv=pinv(Casym); %Gettign the pseudoinverse of the EMGreactive and EMGcontext
Wasym = Cinv*Ymodel; %x= y/C
hatEMGasym_Pos=  Casym * Wasym  ; %yhat = C 

%LMSE using negative pertubation
reactive=find(strcmp(epochOfInterest,'Ramp')==1);
context= find(strcmp(epochOfInterest,'Optimal')==1);
Casym=[C(:,reactive) C(:,context)]; 
Cinv=pinv(Casym); %Gettign the pseudoinverse of the EMGreactive and EMGcontext
Wasym = Cinv*Ymodel; %x= y/C
hatEMGasym_PosP=  Casym * Wasym  ; %yhat = C 


% LMSE C using adaptation PCA
YA=Yasym(35:480,:);
[pp_A,cc_A,aa_A]=pca(YA);%,'Centered','off');
%%Input has to be row observation  and columns variables 
%%pp - Corresponding matrix of eigenvectors
%%cc - data projected in the principal component 
%%aa - vector of eigent values
% Ymean=nanmean(Yasym,1);
Casym = pp_A(:,1:2);
Cinv=pinv(Casym); %Gettign the pseudoinverse of the EMGreactive and EMGcontext
Wasym_A = Cinv*Ymodel; %x= y/C
hatEMGasym_A=  Casym * Wasym_A  ; %yhat = C 


% LMSE C using 5 prior +30 next stage  PCA
YA=[Yasym(1:70,:); Yasym(450:481+60,:)];


% NNMF
swift=abs(min(Yasym',[],'all'));
data=YA+swift;
data2=Yasym+swift;
[W,H] = nnmf(data,4);
data_hat=W*H;
data_hat2=data_hat-swift;

Casym = H';
Cinv=pinv(Casym); %Gettign the pseudoinverse of the EMGreactive and EMGcontext
Wasym_NNMF = Cinv*data2'; %x= y/C
Wasym_NNMF_noshift = Cinv*Ymodel; %x= y/C
hatEMGasym_NNMF=  (Casym * Wasym_NNMF) - swift ; %yhat = C 
hatEMGasym_NNMF_noshift=  (Casym * Wasym_NNMF_noshift); %yhat = C 


[pp_A,cc_A,aa_A]=pca(YA,'Centered','off');
%%Input has to be row observation  and columns variables 
    %%pp - Corresponding matrix of eigenvectors
%%cc - data projected in the principal component 
%%aa - vector of eigent values
% Ymean=nanmean(Yasym,1);
Casym = pp_A(:,1:2);
Cinv=pinv(Casym); %Gettign the pseudoinverse of the EMGreactive and EMGcontext
dYasm=Ymodel-mean(Ymodel,2);
Wasym_A = Cinv*Ymodel; %x= y/C
hatEMGasym_S=  (Casym * Wasym_A); %yhat = C 
 

% LMSE C using post-adaptation PCA
YP=Yasym(481:end,:);
[pp_P,cc_P,aa_P]=pca(YP);%,'Centered','off');
%%Input has to be row observation  and columns variables 
% Ymean=nanmean(Ypost,1);
Casym = pp_P(:,1:2);
Cinv=pinv(Casym); %Gettign the pseudoinverse of the EMGreactive and EMGcontext
Wasym_P = Cinv*Ymodel; %x= y/C
hatEMGasym_P=  Casym * Wasym_P  ; %yhat = C 


%C from PCA 
[pp_2,cc_,aa_]=pca(Ymodel');%,'Centered','off');
%%Input has to be row observation  and columns variables 
Ymean=nanmean(Yasym,1);
Xcentered = cc_(:,1:2)*pp_2(:,1:2)'; 
Ynew=Xcentered+Ymean;


%C from all data LMSE
Casym = pp_2(:,1:2);
Cinv=pinv(Casym); %Gettign the pseudoinverse of the EMGreactive and EMGcontext
Wasym_P = Cinv*Ymodel; %x= y/C
hatEMGasym_all=  Casym * Wasym_P  ; %yhat = C 

% subplot(3,1,3)
%%
binw=5;
figure 
hold on
ea=681
la=880

aux2= 1 - sum((Yasym(ea:la,:)'-hatEMGasym(:,ea:la)).^2)./sum((Yasym(ea:la,:)'- mean(Yasym(ea:la,:)')).^2);
aux2=conv(aux2,ones(1,binw)/binw,'valid'); %Smoothing
plot(aux2,'LineWidth',2,'DisplayName','Linear regression - Ramp + Optimal + Removal','Color','m') ;

aux2= 1 - sum((Yasym(ea:la,:)'-hatEMGasym_S(:,ea:la)).^2)./sum((Yasym(ea:la,:)'- mean(Yasym(ea:la,:)')).^2);
aux2=conv(aux2,ones(1,binw)/binw,'valid'); %Smoothing
% plot(aux2,'LineWidth',2,'DisplayName','PCA+LSE','Color',"#0072BD") ;


aux2= 1 - sum((Yasym(ea:la,:)'-hatEMGasym_Neg(:,ea:la)).^2)./sum((Yasym(ea:la,:)'- mean(Yasym(ea:la,:)')).^2);
aux2=conv(aux2,ones(1,binw)/binw,'valid'); %Smoothing
plot(aux2,'LineWidth',2,'DisplayName','Linear Regression - negative + optimal','Color',"#0072BD") ;

aux2= 1 - sum((Yasym(ea:la,:)'-hatEMGasym_Neg2(:,ea:la)).^2)./sum((Yasym(ea:la,:)'- mean(Yasym(ea:la,:)')).^2);
aux2=conv(aux2,ones(1,binw)/binw,'valid'); %Smoothing
plot(aux2,'LineWidth',2,'DisplayName','Linear Regression - negative + optimal + ramp','Color',"#77AC30") ;


aux2= 1 - sum((Yasym(ea:la,:)'-hatEMGasym_Pos(:,ea:la)).^2)./sum((Yasym(ea:la,:)'- mean(Yasym(ea:la,:)')).^2);
aux2=conv(aux2,ones(1,binw)/binw,'valid'); %Smoothing
plot(aux2,'LineWidth',2,'DisplayName','Linear Regression - Removal + Optimal','Color',"g") ;

aux2= 1 - sum((Yasym(ea:la,:)'-hatEMGasym_PosP(:,ea:la)).^2)./sum((Yasym(ea:la,:)'- mean(Yasym(ea:la,:)')).^2);
aux2=conv(aux2,ones(1,binw)/binw,'valid'); %Smoothing
plot(aux2,'LineWidth',2,'DisplayName','Linear Regression -Ramp + Optimal','Color',"k") ;

% aux2= 1 - sum((Yasym'-hatEMGasym_A).^2)./sum((Yasym'- mean(Yasym')).^2);
% aux2=conv(aux2,ones(1,binw)/binw,'valid'); %Smoothing
% plot(aux2,'LineWidth',2,'DisplayName','PCA_{base+adaptation}','Color',"#0072BD") ;
%  
% 
% aux2= 1 - sum((Yasym'-hatEMGasym_P).^2)./sum((Yasym'- mean(Yasym')).^2);
% aux2=conv(aux2,ones(1,binw)/binw,'valid'); %Smoothing
% plot(aux2,'LineWidth',2,'DisplayName','PCA_{Post}','Color',"#77AC30") ;

aux2= 1 - sum((Yasym'-Ynew').^2)./sum((Yasym'- mean(Yasym')).^2);
aux2=conv(aux2,ones(1,binw)/binw,'valid'); %Smoothing
% plot(aux2,'LineWidth',2,'DisplayName','PCA_{all}','Color',"#4DBEEE") ;


aux2= 1 - sum((Yasym'-hatEMGasym_all).^2)./sum((Yasym'- mean(Yasym')).^2);
aux2=conv(aux2,ones(1,binw)/binw,'valid'); %Smoothing
% plot(aux2,'LineWidth',2,'DisplayName','PCA_{LSE}','Color',"#A2142F") ;


% aux2= 1 - sum((Yasym'-hatEMGasym_fit).^2)./sum((Yasym'- mean(Yasym')).^2);
% aux2=conv(aux2,ones(1,binw)/binw,'valid'); %Smoothing
% plot(aux2,'LineWidth',2,'DisplayName','Data_{adapt}','Color',"#EDB120") 
% 
% aux2= 1 - sum((Yasym'-hatEMGasym_Pfit).^2)./sum((Yasym'- mean(Yasym')).^2);
% aux2=conv(aux2,ones(1,binw)/binw,'valid'); %Smoothing
% plot(aux2,'LineWidth',2,'DisplayName','Data_{post}','Color',[255/255 0/255 255/255]) 


% aux2= 1 - sum((data'-data_hat').^2)./sum((data'- mean(data')).^2);
% aux2=conv(aux2,ones(1,binw)/binw,'valid'); %Smoothing
% plot(aux2,'LineWidth',2,'DisplayName','NNMF+swift','Color',[255/255 0/255 255/255]) 



aux2= 1 - sum((Yasym(ea:la,:)'-hatEMGasym_NNMF(:,ea:la)).^2)./sum((Yasym(ea:la,:)'- mean(Yasym(ea:la,:)')).^2);
aux2=conv(aux2,ones(1,binw)/binw,'valid'); %Smoothing
plot(aux2,'LineWidth',2,'DisplayName','NNMF+LSE','Color','r') 

aux2= 1 - sum((Yasym'-hatEMGasym_NNMF_noshift).^2)./sum((Yasym'- mean(Yasym')).^2);
aux2=conv(aux2,ones(1,binw)/binw,'valid'); %Smoothing
% plot(aux2,'LineWidth',2,'DisplayName','NNMF+LSE (not shifted data)','Color','b') 

% aux2= 1 - sum((Yasym'-data_hat2').^2)./sum((Yasym'- mean(Yasym')).^2);
% aux2=conv(aux2,ones(1,binw)/binw,'valid'); %Smoothing
% plot(aux2,'LineWidth',2,'DisplayName','NNMF-swift','Color','g') 



% legend('AutoUpdate','off') 
ylabel('R^2')
  pp=patch([0 440  440 0],[0 0 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
 
% Xcentered = cc_(:,1:3)*pp_2(:,1:3)'; 
% % Ymean=nanmean(Yasym(41:480,:),1);
% Ynew=Xcentered+Ymean;
% aux2= 1 - sum((Yasym'-Ynew').^2)./sum((Yasym'- mean(Yasym')).^2);
% aux2=conv(aux2,ones(1,binw)/binw,'valid'); %Smoothing
% plot(aux2,'LineWidth',2,'DisplayName','PCA_{3}','Color','r') ;
legend('AutoUpdate','off') 
ylabel('R^2')
%% Regression
GreenOrangeGradient
 yt=1:14;
 fs=14;
ytl={'GLU','TFL','HIP','RF','VL','VM','BF', 'SEMB','SEMT','MG','LG','SOL','PER','TA'};
% set(gca,'XTick',[],'YTick',yt,'YTickLabel',ytl,'FontSize',fs)


figure
subplot(2,2,1)
imagesc((reshape(C_dulce(:,1),12,14)'),[-.5 .5]);
colormap(flipud(map))
colorbar
set(gca,'XTick',[],'YTick',yt,'YTickLabel',ytl,'FontSize',fs)
title('Reactive')

subplot(2,2,2)
plot(W_dulce(1,:))
axis tight
xlabel('strides')
ylabel('W')
 
subplot(2,2,3) 
imagesc((reshape(C_dulce(:,2),12,14)'),[-.5 .5]);
colormap(flipud(map))
colorbar
set(gca,'XTick',[],'YTick',yt,'YTickLabel',ytl,'FontSize',fs)
title('Context-Specific')

subplot(2,2,4)
plot(W_dulce(2,:))
axis tight  
xlabel('strides')
ylabel('W')
set(gcf,'color','w');

%% NNMF
GreenOrangeGradient
 yt=1:14;
ytl={'GLU','TFL','HIP','RF','VL','VM','BF', 'SEMB','SEMT','MG','LG','SOL','PER','TA'};
% set(gca,'XTick',[],'YTick',yt,'YTickLabel',ytl,'FontSize',fs)

figure
subplot(2,2,1)
ht=H';
imagesc((reshape(ht(:,1),12,14)'),[0 0.16]);
colormap(flipud(map))
colorbar
set(gca,'XTick',[],'YTick',yt,'YTickLabel',ytl,'FontSize',fs)
title('Factor 1')

subplot(2,2,2)
plot(Wasym_NNMF(1,:))
axis tight
xlabel('strides')
ylabel('W')
 
subplot(2,2,3)
imagesc((reshape(ht(:,2),12,14)'),[0 0.16]);
colormap(flipud(map))
colorbar
set(gca,'XTick',[],'YTick',yt,'YTickLabel',ytl,'FontSize',fs)
title('Factor 2')

subplot(2,2,4)
plot(Wasym_NNMF(2,:))
axis tight  
xlabel('strides')
ylabel('W')

set(gcf,'color','w');
%% NNMF no shited 
GreenOrangeGradient
 yt=1:14;
ytl={'GLU','TFL','HIP','RF','VL','VM','BF', 'SEMB','SEMT','MG','LG','SOL','PER','TA'};
% set(gca,'XTick',[],'YTick',yt,'YTickLabel',ytl,'FontSize',fs)

figure
subplot(2,2,1)
ht=H';
imagesc((reshape(ht(:,1),12,14)'),[0 0.16]);
colormap(flipud(map))
colorbar
set(gca,'XTick',[],'YTick',yt,'YTickLabel',ytl,'FontSize',fs)
title('Factor 1')

subplot(2,2,2)
plot(Wasym_NNMF_noshift(1,:))
axis tight
xlabel('strides')
ylabel('W')
title('Factor 1 - Data no shifted ')
     
subplot(2,2,3)
imagesc((reshape(ht(:,2),12,14)'),[0 0.16]);
colormap(flipud(map))
colorbar
set(gca,'XTick',[],'YTick',yt,'YTickLabel',ytl,'FontSize',fs)
title('Factor 2')

subplot(2,2,4)
plot(Wasym_NNMF_noshift(2,:))
axis tight  
xlabel('strides')
ylabel('W')
title('Factor 2 - Data no shifted ')

set(gcf,'color','w');

%%

mm= 0:12:168;
mm2=1:12:168;
for m=1:14
    
    muscles_regression(m,:)= 1 - sum((Ymodel(mm2(m):mm(m+1),:)-hatEMGasym(mm2(m):mm(m+1),:)).^2)./sum((Ymodel(mm2(m):mm(m+1),:)- mean(Ymodel(mm2(m):mm(m+1),:))).^2);
    muscles_nnmf(m,:)= 1 - sum((Ymodel(mm2(m):mm(m+1),:)-hatEMGasym_NNMF(mm2(m):mm(m+1),:)).^2)./sum((Ymodel(mm2(m):mm(m+1),:)- mean(Ymodel(mm2(m):mm(m+1),:))).^2);
 
end


figure
title('R^2')
ytl={'GLU','TFL','HIP','RF','VL','VM','BF', 'SEMB','SEMT','MG','LG','SOL','PER','TA'};
for m=1:7
    subplot(7,1,m)
    hold on
    plot(movmean(muscles_regression(m,:),5))
    plot(movmean(muscles_nnmf(m,:),5))
%     scatter(1:length(muscles(m,:)),movmean(muscles(m,:),5),'filled')
    legend({'regression','NNMF'},'AutoUpdate','off');
    ylabel(ytl{m})
     yline(0)
    
%     pp=patch([40 480 480 40],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
   
    pp=patch([41 480  480 41],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
 ylim([-1 1])
uistack(pp,'bottom')
end
set(gcf,'color','w')

figure
title('R^2')
for m=8:14
    subplot(7,1,m-7)
    hold on
        plot(movmean(muscles_regression(m,:),5))
        plot(movmean(muscles_nnmf(m,:),5))
%     scatter(1:length(muscles(m,:)),movmean(muscles(m,:),5),'filled')
    legend(ytl{m},'AutoUpdate','off');
    legend({'regression','NNMF'},'AutoUpdate','off');
    ylabel(ytl{m})
    yline(0)
%     yline(nanmean(muscles(m,10:30)))
%     pp=patch([40 480 480 40],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
    pp=patch([41 480  480 41],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
    ylim([-1 1])
uistack(pp,'bottom')
end
set(gcf,'color','w')


%% ONLY adaptation PCA

YA=Yasym(41:480,:);
[pp_A,cc_A,aa_A]=pca(YA);%,'Centered','off');
%%Input has to be row observation  and columns variables 
%%pp - Corresponding matrix of eigenvectors
%%cc - data projected in the principal component 
%%aa - vector of eigent values
% Ymean=nanmean(Yasym,1);
Ymean=nanmean(YA,1);
Xcentered = cc_A(:,1:2)*pp_A(:,1:2)'; 
Ynew=Xcentered+Ymean;


binw=5;
figure
aux1= 1 - sum((YA'-Ynew').^2)./sum((YA'- mean(YA')).^2);
aux1=conv(aux1,ones(1,binw)/binw,'valid'); %Smoothing
plot(aux1,'LineWidth',2,'DisplayName','PCA_{adapt}','Color',[255/255 0/255 255/255]) 

%% ONLY post-adaptation PCA

YP=Yasym(481:end,:);
[pp_P,cc_P,aa_P]=pca(YP);%,'Centered','off');
%%Input has to be row observation  and columns variables 
Ymean=nanmean(YP,1);
Xcentered = cc_P(:,1:2)*pp_P(:,1:2)'; 
Ynew=Xcentered+Ymean;

binw=5;
figure
aux2= 1 - sum((YP'-Ynew').^2)./sum((YP'- mean(YP')).^2);
aux2=conv(aux2,ones(1,binw)/binw,'valid'); %Smoothing
plot(aux2,'LineWidth',2,'DisplayName','PCA_{post}','Color',[255/255 0/255 255/255]) 



legend('AutoUpdate','off') 
ylabel('R^2')
%%
% figure
aux3=[nan(1,40) aux1 aux2]
plot(aux3,'LineWidth',2,'DisplayName','PCA_concatenated','Color','r') 

%   pp=patch([40 480  480 40],[-.2 -.2 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');

% [pp_3,cc_3,aa3]=pca(Yasym);
%%
%Recontruction of the data 
Xcentered = cc_(:,1:2)*pp_2(:,1:2)'; 

Ynew=Xcentered+Ymean;


figure 
subplot(3,1,1)
plot(aa_,'LineWidth',2)
ylabel('Eigenvalues')
xline(2,'r')
xlabel('Number of Components')
% ylim([-.5 7])
 axis tight

% yline(0.85,'r')
% xline(0.85,'r')

subplot(3,1,2)
csum=cumsum(aa_);
variance_explained = csum / sum(aa_);
plot(variance_explained,'LineWidth',2)
xline(2,'r')
ylabel('Variance explained')
xlabel('Number of Components')
set(gcf,'color','w')


binw=5;
subplot(3,1,3)
hold on 
aux2= 1 - sum((Yasym'-Ynew').^2)./sum((Yasym'- mean(Yasym')).^2);
aux2=conv(aux2,ones(1,binw)/binw,'valid'); %Smoothing
plot(aux2,'LineWidth',2,'DisplayName','PCA_{2}','Color','k') ;


Xcentered = cc_(:,1:3)*pp_2(:,1:3)'; 
% Ymean=nanmean(Yasym(41:480,:),1);
Ynew=Xcentered+Ymean;
aux2= 1 - sum((Yasym'-Ynew').^2)./sum((Yasym'- mean(Yasym')).^2);
aux2=conv(aux2,ones(1,binw)/binw,'valid'); %Smoothing
plot(aux2,'LineWidth',2,'DisplayName','PCA_{3}','Color','r') ;
%% Using first 2 components to estimate the states 

C=[pp_2(:,1:3)];
Cinv=pinv(C);

X= Yasym*Cinv';

figure
subplot(3,1,1)
plot(movmean(X(:,1),5),'LineWidth',2)
yline(0)
title('C_1 = PC1')
% legend('ATR','ATS')

subplot(3,1,2)
plot(movmean(X(:,2),5),'LineWidth',2)
yline(0)
title('C_2 = PC2')
% 
    subplot(3,1,3)
    plot(movmean(X(:,3),5),'LineWidth',2)
    yline(0)
    title('C_3 = PCA3')
    legend('ATR','ATS')

set(gcf,'color','w')
%%
model.C=C;
legacy_vizSingleModel_FreeModel_ShortAdaptation(model,Yasym',Uf,0)
% legacy_vizSingleModelMLMC_FreeModel(model,datSet.out,datSet.in)
% 

% model.C=C(:,1:2);
% legacy_vizSingleModelMLMC_FreeModel(model,datSet.out,datSet.in)