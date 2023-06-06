% fname='dynamicsData_BATR_subj_12_RemoveBadMuscles1_splits_0_V4.h5';
fname='dynamicsData_BATS_subj_12_RemoveBadMuscles1_splits_0_WithPost2V2.h5'
% load BATR_12_AsymC16_ShortPertubations_RemovedBadMuscle_1.mat

% ALL 24 participatns
% load BAT_24_AsymC16_ShortPertubations_RemovedBadMuscle_1.mat
load BAT_24_IndvLegsC16_ShortPertubations_RemovedBadMuscle_1.mat


EMGdata=h5read(fname,'/EMGdata');
 
binwith=5;
[Y,Yasym,~,U,~,Ysum,Yinv]=groupDataToMatrixForm_Update(1:size(EMGdata,3),fname,0);
Yasym=Yinv;

Ymodel=Yasym';
swift=abs(min(Yasym',[],'all'));
%LMSE Dulce's fit 
reactive=find(strcmp(epochOfInterest,'Ramp')==1);
reactive2=find(strcmp(epochOfInterest,'Tied post ramp')==1);
context= find(strcmp(epochOfInterest,'Optimal')==1);
base= find(strcmp(epochOfInterest,'TM base')==1);

NNMF=0
PCA_analysis=0
removemean=0
const=0
removebaseline=0

Casym=[C(:,reactive) C(:,context) C(:,reactive2) C(:,base)]; % EMGreactive and EMGcontext
% Casym=[ C(:,base)]; % EMGreactive and EMGcontext

Ymodel=Yasym';
% swift=abs(min(Yasym',[],'all'));
% %LMSE Dulce's fit 
% 
% reactive=find(strcmp(epochOfInterest,'Ramp')==1);
% context= find(strcmp(epochOfInterest,'Optimal')==1);
% reactive2=find(strcmp(epochOfInterest,'Tied post ramp')==1);
% % Casym=[C(:,5) C(:,6)];%./vecnorm([C(:,5) C(:,6)]); 
% Casym=[C(:,reactive) C(:,context) C(:,reactive2) ]; 
% Casym=[C(:,reactive)]; 
Cinv=pinv(Casym); %Gettign the pseudoinverse of the EMGreactive and EMGcontext
Wasym = Cinv*(Ymodel); %x= y/C
hatEMGasym=  (Casym * Wasym); %yhat = C 

C_dulce=Casym;
W_dulce=Wasym;


% LMSE C using 5 prior +30 next stage  PCA
YA=[Yasym(1:70,:); Yasym(450:481+60,:)];


% NNMF
swift=abs(min(Yasym',[],'all'));
data=YA+swift;
data2=Yasym+swift;
[W,H] = nnmf(data,3);
data_hat=W*H;
data_hat2=data_hat-swift;

Casym = H';
Cinv=pinv(Casym); %Gettign the pseudoinverse of the EMGreactive and EMGcontext
Wasym_NNMF = Cinv*data2'; %x= y/C
Wasym_NNMF_noshift = Cinv*Ymodel; %x= y/C
hatEMGasym_NNMF=  (Casym * Wasym_NNMF) - swift ; %yhat = C 
hatEMGasym_NNMF_noshift=  (Casym * Wasym_NNMF_noshift); %yhat = C 

%%
binw=5;
% figure 

ea=1
la=680
figure
aux2= 1 - sum((Yasym(ea:la,:)'-hatEMGasym(:,ea:la)).^2)./sum((Yasym(ea:la,:)'- mean(Yasym(ea:la,:)')).^2);
aux2=conv(aux2,ones(1,binw)/binw,'valid'); %Smoothing
plot(aux2,'LineWidth',2,'DisplayName','Linear regression') ;

hold on

aux2= 1 - sum((Yasym(ea:la,:)'-hatEMGasym_NNMF(:,ea:la)).^2)./sum((Yasym(ea:la,:)'- mean(Yasym(ea:la,:)')).^2);
aux2=conv(aux2,ones(1,binw)/binw,'valid'); %Smoothing
plot(aux2,'LineWidth',2,'DisplayName','NNMF','Color','r') 

ylabel('R^2')
xlabel('Strides')

pp=patch([40 480 480 40],[0.8 0.8 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
 legend({'Linear Regression';'NNMF'},'AutoUpdate','off')
