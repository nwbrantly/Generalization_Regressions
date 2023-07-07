%model fit per muscle

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
clear all;clc
%% Free model - Linear regression - Asymmetry with baseline

% %Testing 
% fname='dynamicsData_BATR_subj_12_RemoveBadMuscles1_splits_0_V4.h5';
% fname='dynamicsData_BAT_subj_24_RemoveBadMuscles1_splits_0_V4.h5'
% fname='dynamicsData_BATR_subj_12_RemoveBadMuscles1_splits_0_WithPost2V2.h5'
% fname='dynamicsData_BATS_subj_12_RemoveBadMuscles1_splits_0_WithPost2V2_WogBaseline.h5'

%posterior muscles
fname='dynamicsData_BATR_subj_12_RemoveBadMuscles1_splits_0_PosteriorMuscles.h5';
load BAT_24_IndvLegsC17_ShortPertubations_RemovedBadMuscle_1RemoveBias_0_PosteriorMuscles.mat

%BATR 
% load BATR_12_AsymC16_ShortPertubations_RemovedBadMuscle_1.mat
% load BATR_12_IndvLegsC16_ShortPertubations_RemovedBadMuscle_1RemovBias_0.mat

%BATS
% load BATS_12_IndvLegsC17_ShortPertubations_RemovedBadMuscle_1RemovBias_0.mat

% ALL 24 participatns
% load BAT_24_AsymC16_ShortPertubations_RemovedBadMuscle_1.mat
% load BAT_24_IndvLegsC16_ShortPertubations_RemovedBadMuscle_1.mat


EMGdata=h5read(fname,'/EMGdata');
 
binwith=5;
[Y,Yasym,~,U,~,Ysum,Yinv,labels]=groupDataToMatrixForm_Update(1:size(EMGdata,3),fname,0);
Yasym=Yinv;
% Yasym=Yasym(:,1:size(Yasym,2)/2); %SLOW
% Yasym=Yasym(:,169:end); %FAST

% range=61:72;
% Yasym=Yasym(:,range); %One mucles
% Yasym=Yinv(1:480,:);
% Yasym=Y;
% [Y,Yasym,Ycom,U,Ubreaks,Ysum]
% Yasym=Yasym(481:end,:);
% U=U(481:end);
Uf=[U;ones(size(U))];


% C=[C1 C3];
% C=Cnew;
reactive=find(strcmp(epochOfInterest,'Ramp')==1);
reactive2=find(strcmp(epochOfInterest,'Tied post ramp')==1);
context= find(strcmp(epochOfInterest,'Optimal')==1);
base= find(strcmp(epochOfInterest,'TM base')==1);
% base= find(strcmp(epochOfInterest,'OG base')==1);
% Casym=[C(:,reactive) C(:,context) -C(:,reactive2) C(:,base)]; % EMGreactive and EMGcontext

reactive2=find(strcmp(epochOfInterest,'NegShort_{late}')==1);
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
    Ymodel=Ymodel-bias'; % Transpose the EMG data to match equations
%     Ymodel=Ymodel(:,481:680);

        Cmuscles=reshape(Casym',2,12,12);
        Ymuscles=reshape(Ymodel',880,12,12);
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


    
    %% STEPWISE regression
    clear temp temp2 temp3
    for i=1:12
        %reactive 
        temp=pinv(Cmuscles(1,:,i)');
        Xreactive =temp*Ymuscles(:,:,i)'; %x= y/C
        Yestimated=  Cmuscles(1,:,i)'* Xreactive ;
        residual=Ymuscles(:,:,i)-Yestimated';
        
        %contextual
        temp3=pinv(Cmuscles(2,:,i)');
        Xcontextual =temp3*residual'; %x= y/C
        Yestimated2=  Cmuscles(2,:,i)'* Xcontextual;
        residual2=residual-Yestimated2';
        
        
        temp2(:,:,i)=[Xreactive'  Xcontextual'];
        model{i}.C=Cmuscles(:,:,i);
        
        
    end

    
end
Uf=Uf(:,1:size(temp2,1));

% t=12:-1:1;
%%
analysis=0
for i=4
    Xasym=[temp2(1:40,:,i);nan(1,size(temp2,2));temp2(41:480,:,i);nan(1,size(temp2,2));temp2(481:end,:,i)];
%      Xasym=[temp2(1:200,:,i)];
    % Xasym=[X2asym(1:40,:);nan(1,size(X2asym,2));X2asym(41:80,:);...
    %     nan(1,size(X2asym,2));X2asym(81:520,:);nan(1,size(X2asym,2));X2asym(521:end,:)];
    % Xasym=[X2asym(681:843,:);nan(1,size(X2asym,2))];
    
    figure
    subplot(2,1,1)
    hold on
    % scatter(1:length(movmean(Xasym(:,1),binwith)), movmean(Xasym(:,1),binwith),10,'k','filled')
    scatter(1:length(movmean(Xasym(:,1),binwith)), movmean(Xasym(:,1),binwith),'filled','MarkerFaceColor',"#EDB120") %"#77AC30" )%
%     pp=patch([40 480 480 40],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
    % plot( movmean(Xasym(:,1),binwith),'Color',"#EDB120",'LineWidth',5)
    % pp=patch([40 480 480 40],[0.5 0.5 1.5 1.5],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
    % pp=patch([40 480 480 40],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
    % pp=patch([0 440 440 0],[-0.5 -0.5 1 1],.7*ones(1reactive2,3),'FaceAlpha',.2,'EdgeColor','none');
    
    % pp=patch([40 480 480 40],[-4 -4 5 5],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
    
    % legend('Baseline','AutoUpdate','off')
    legend('Reactive','AutoUpdate','off')
    title(labels(i).Data)
    %     legend('Removal Perturbation','AutoUpdate','off')
    % uistack(pp,'bottom')
    yline(0)
    ylabel({'Reactive';'(A.U)'})
    %     ylabel({'Removal';'(A.U)'})
    xlabel('strides')
    
    if size(temp2(:,:,i),2)>=2
        % figure
        subplot(2,1,2)
        hold on
        scatter(1:length(movmean(Xasym(:,2),binwith)), movmean(Xasym(:,2),binwith),'filled','MarkerFaceColor',"#77AC30")
%         pp=patch([40 480 480 40],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
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
    if i<=6
    isF=0;
    else
      isF=1;  
    end
    legacy_vizSingleModel_FreeModel_ShortAdaptation_IndvStepWise(model{i},Ymuscles(:,:,i)',Uf,analysis,{labels(i).Data(2:end-1)},isF)
    set(gcf,'color','w')
end

%% Define the type of analysis that you want to do 
%0 linear regression 
%1 Using hypothesize patterns as constants
%2 PCA 
%3 removing the mean for the data (varaince reconstrucitons)

    
 

% legacy_vizSingleModel_FreeModel_ShortAdaptation(model,Ymodel,Uf,analysis)


% legacy_vizSingxleModel_FreeModel_›ShortAdae ptation(model,Ymodel',Uf,1)
%     
% model.C=C(:,1); 
% legacy_vizSingleModel_FreeeModel_ShortAdaptation(model,Ymodel',Uf,1)

% legacy_vizSingleModel_FreeModel_ShortAdaptation_IndvLeg_cond(model,Ymodel,Uf,analysis)
