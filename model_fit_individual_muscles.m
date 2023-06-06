%Model fit per muscle

%% Adding path 

addpath(genpath('/Users/dulcemariscal/Documents/GitHub/Generalization_Regressions'))
addpath(genpath('/Users/dulcemariscal/Documents/GitHub/labTools'))
addpath(genpath('/Users/dulcemariscal/Documents/GitHub/LongAdaptation'))
addpath(genpath('/Users/dulcemariscal/Documents/GitHub/R01'))
addpath(genpath('/Users/dulcemzariscal/Documents/GitHub/splitbelt-EMG-adaptation'))
addpath(genpath('/Users/dulcemariscal/Documents/GitHub/EMG-LTI-SSM'))
addpath(genpath('/Users/dulcemariscal/Documents/GitHub/matlab-linsys'))
rmpath(genpath('/Users/dulcezmariscal/Documents/GitHub/PittSMLlab'))

%% Load real data:
clear all;clc;close all
%% Free model - Linear regression - Asymmetry with baseline
%% This is just the saved data - Update accrodingly 

type= 1; %1 testing; 2 training ; 3 both groups
%EMG Data 
if type==1 %Testing
    fname='dynamicsData_BATS_subj_12_RemoveBadMuscles1_splits_0_WithPost2V2_WogBaseline.h5'
elseif type==2 %Training
    %posterior muscles
    % fname='dynamicsData_BATR_subj_12_RemoveBadMuscles1_splits_0_PosteriorMuscles.h5';
    
    % fname='dynamicsData_BATR_subj_12_RemoveBadMuscles1_splits_0_V4.h5'
    fname='dynamicsData_BATR_subj_12_RemoveBadMuscles1_splits_0_WithPost2V2.h5'
elseif type==3 %Both groups
    fname='dynamicsData_BAT_subj_24_RemoveBadMuscles1_splits_0_V4.h5'
end

%%%%%%%% Regressors

%BATR
% load BATR_12_AsymC16_ShortPertubations_RemovedBadMuscle_1.mat
% load BATR_12_IndvLegsC16_ShortPertubations_RemovedBadMuscle_1RemovBias_0.mat

%BATS
% load BATS_12_IndvLegsC17_ShortPertubations_RemovedBadMuscle_1RemovBias_0.mat

% ALL 24 participants
% load BAT_24_AsymC16_ShortPertubations_RemovedBadMuscle_1.mat
load BAT_24_IndvLegsC16_ShortPertubations_RemovedBadMuscle_1.mat
% load BAT_24_IndvLegsC17_ShortPertubations_RemovedBadMuscle_1RemoveBias_0_PosteriorMuscles.mat

%% Getting the data

EMGdata=h5read(fname,'/EMGdata');
binwith=10;
[Y,Yasym,~,U,~,Ysum,Yinv,labels]=groupDataToMatrixForm_Update(1:size(EMGdata,3),fname,0);
Yasym=Y;
Uf=[U;ones(size(U))];

%% Organizing the data

%Pick the conditions that we are going to use for the regression
reactive=find(strcmp(epochOfInterest,'Ramp')==1);
context= find(strcmp(epochOfInterest,'Optimal')==1);
reactive2=find(strcmp(epochOfInterest,'NegShort_{late}')==1);


removebaseline=1; %Remove bias flag


Casym=[C(:,reactive2) C(:,context)]; % EMGreactive and EMGcontext
Ymodel=Yasym';  % Transpose the EMG data to match equations


if removebaseline==1
    bias=nanmean(Yasym(5:30,:)); %Computing the bias
    
    Casym=Casym-bias'; %Removing bias from the regressors
    Ymodel=Ymodel-bias'; %removing bias from the data
    
    %Organizing the data per muscle
    Cmuscles=reshape(Casym',2,12,size(Yasym,2)/12); %muscles for regressors 2 number of regressors, 12 number of phase of the gait cycle
    Ymuscles=reshape(Ymodel(:,1:680)',680,12,size(Yasym,2)/12); %data
end

%% Linear regression individual muscles
reconstruction_indv=[];
data=[];
C_indv=[];
X_indv=[];

for i=1:28
    unit(:,:,i)=Cmuscles(:,:,i)'./vecnorm(Cmuscles(:,:,i)'); %Getting the unit vector of the regressors
    temp(:,:,i)=pinv(unit(:,:,i)'); % Getting the inverse
    X2asym(:,:,i) =temp(:,:,i)'*Ymuscles(:,:,i)'; %x= y/C
    Y2asym(:,:,i)=  unit(:,:,i)* X2asym(:,:,i) ; %yhat
    temp2(:,:,i)=X2asym(:,:,i)'; % transposing the dynamics vector
    model{i}.C=unit(:,:,i); %saving the regressors. Yhis is necesary for the plotting funciton
    
    reconstruction_indv =[ reconstruction_indv ; Y2asym(:,:,i)]; % Concatenating the data reconstructed
    data =[ data ; Ymuscles(:,:,i)'];  % Concatenating the data
    C_indv=[C_indv;unit(:,:,i)];  % Concatenating the regressors for each muscles
    X_indv=[X_indv,X2asym(:,:,i) ];  % Concatenating the dynamics
    
    % Checking for colinearity and correlation between the regresssions
    temp5=corrcoef(model{i}.C); % correlation coeficient
    correlation(i,1)=temp5(2);
    temp10(i,:)=vif([model{i}.C nanmean(Ymuscles(481:491,:,i))']); %Variance inflation
    temp5=vif([model{i}.C]);
    impact(i,:)=temp5;
    
end

index=find(impact>5); % fidn the regressor with high colinearity

Uf=Uf(:,1:size(temp2,1));
set(gcf,'color','w')
Uf=Uf(:,1:size(temp2,1));

%%  Plotting 
model2{1}.C=C_indv; %Save regressors for model
model2{1}.Out= reconstruction_indv; %Save reconstruction of the indiivudal muscles to the model format for plotting 
analysis=0; %Flag asking if you want to run the analysis
isF=0; %Flag for fast leg


legacy_vizSingleModel_FreeModel_ShortAdaptation_IndvLeg2(model2{1},data,Uf,analysis,[],isF)
% legacy_vizSingleModel_FreeModel_ShortAdaptation_IndvLeg(model2{1},data,Uf,analysis,[],isF)
% save('BATS_indv_muscles.mat','X2asym')


%% Plotting time course for individual muscles
analysis=0

for i=28 %pick the muscle that you want
    
    % Pick the data that you want to plot 
    %    Xasym=[temp2(1:40,:,i);nan(1,size(temp2,2));temp2(41:480,:,i);nan(1,size(temp2,2));temp2(481:end,:,i)];
    Xasym=[temp2(481:end,:,i)];
    %      Xasym=[temp2(1:200,:,i)];
    % Xasym=[X2asym(1:40,:);nan(1,size(X2asym,2));X2asym(41:80,:);...
    %     nan(1,size(X2asym,2));X2asym(81:520,:);nan(1,size(X2asym,2));X2asym(521:end,:)];
    % Xasym=[X2asym(681:843,:);nan(1,size(X2asym,2))];
    
    figure
    subplot(2,1,1)
    hold on
    % scatter(1:length(movmean(Xasym(:,1),binwith)), movmean(Xasym(:,1),binwith),10,'k','filled')
    scatter(1:length(movmean(Xasym(:,1),binwith)), movmean(Xasym(:,1),binwith),'filled','MarkerFaceColor',"#EDB120") %"#77AC30" )%
        
    % legend('Baseline','AutoUpdate','off')
    legend('Negative')
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
        scatter(1:length(movmean(Xasym(:,2),binwith)), movmean(Xasym(:,2),binwith),'filled','MarkerFaceColor'," #00008B")
        %
        legend('Contextual')
        % % legend('Switch','AutoUpdate','off')
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
    % Plot the time course plus the regressors 
    legacy_vizSingleModel_FreeModel_ShortAdaptation_IndvLeg(model{i},Ymuscles(:,:,i)',Uf,analysis,{labels(i).Data(2:end-1)},isF)
    set(gcf,'color','w')
    
    
end

%% 
%% Geting the average of the first 5 strides post-adaptation per muscle
figure
hold on
for i=1:28
    if i<15
        Li{1}=scatter(nanmean(temp2(481:485,1,i)),nanmean(temp2(481:485,2,i)),100,"filled",'MarkerFaceColor', 'b');
        text(nanmean(temp2(481:485,1,i)),nanmean(temp2(481:485,2,i)),{labels(i).Data(2:end-1)})
    else
        Li{2}=scatter(nanmean(temp2(481:485,1,i)),nanmean(temp2(481:485,2,i)),100,"filled",'MarkerFaceColor', 'r')  ;
        text(nanmean(temp2(481:485,1,i)),nanmean(temp2(481:485,2,i)),{labels(i).Data(2:end-1)})
    end
    
end
legend([Li{:}],['Slow';'Fast'])
ylabel({'Contextual';'A.U'})
xlabel({'Reactive';'A.U'})

xlim([-0.3 1.4])
ylim([-0.3 1.4])
set(gcf,'color','w')

%% Geeting the average of Early Adapt and lateAdapt
figure
hold on
for i=1:28
    if i<15
        Li{1}=scatter(nanmean(temp2(41:45,1,i)),nanmean(temp2(41:45,2,i)),100,"filled",'MarkerFaceColor', 'b');
        text(nanmean(temp2(41:45,1,i)),nanmean(temp2(41:45,2,i)),{labels(i).Data(2:end-1)})
    else
        Li{2}=scatter(nanmean(temp2(41:45,1,i)),nanmean(temp2(41:45,2,i)),100,"filled",'MarkerFaceColor', 'r')  ;
        text(nanmean(temp2(41:45,1,i)),nanmean(temp2(41:45,2,i)),{labels(i).Data(2:end-1)})
    end
    
end
legend([Li{:}],['Slow';'Fast'])
ylabel({'Contextual';'A.U'})
xlabel({'Reactive';'A.U'})
% axis square
xlim([-1 3])
ylim([-1 1])

title('EarlyAdapt')
set(gcf,'color','w')


figure
hold on
for i=1:28
    if i<15
        Li{1}=scatter(nanmean(temp2(440:480,1,i)),nanmean(temp2(440:480,2,i)),100,"filled",'MarkerFaceColor', 'b');
        text(nanmean(temp2(440:480,1,i)),nanmean(temp2(440:480,2,i)),{labels(i).Data(2:end-1)})
    else
        Li{2}=scatter(nanmean(temp2(440:480,1,i)),nanmean(temp2(440:480,2,i)),100,"filled",'MarkerFaceColor', 'r')  ;
        text(nanmean(temp2(440:480,1,i)),nanmean(temp2(440:480,2,i)),{labels(i).Data(2:end-1)})
    end
    
end
legend([Li{:}],['Slow';'Fast'])
ylabel({'Contextual';'A.U'})
xlabel({'Reactive';'A.U'})
axis square
title('LateAdapt')
set(gcf,'color','w')
xlim([-.1 .85])
ylim([-.1 .85])

%% Comparing individual analysis with the leg specific analysis

Cslow=[Casym(1:size(C,1)/2,1)  Casym(1:size(C,1)/2,2)]./vecnorm([Casym(1:size(C,1)/2,1)  Casym(1:size(C,1)/2,2)]); %SLOW
Cfast=[Casym(1+size(C,1)/2:end,1)  Casym(1+size(C,1)/2:end,2)]./vecnorm([Casym(1+size(C,1)/2:end,1)  Casym(1+size(C,1)/2:end,2)]); %FAST

Ymodel=Ymodel';
Yslow=Ymodel(1:680,1:size(Ymodel,2)/2); %SLOW
Yfast=Ymodel(1:680,size(Ymodel,2)/2+1:end); %FAST



Cinv=pinv(Cslow);
Xslow = Cinv*Yslow'; %x= y/C
hatYslow=  Cslow *Xslow ; %yhat = C

% hatYslow= hatYslow';

Cinv=pinv(Cfast);
Xfast = Cinv*Yfast'; %x= y/C
hatYfast=  Cfast *Xfast ; %yhat = C

% hatYfast=hatYfast;
Yfast=Yfast';
Yslow=Yslow';

% Variance explained

ex2=[0.2314    0.2980    0.7529];
ex1=[0.7255    0.0863    0.1608];
mid=ones(1,3);
N=100;
gamma=1.5; %gamma > 1 expands the white (mid) part of the map, 'hiding' low values. Gamma<1 does the opposite
gamma=1;
map=[flipud(mid+ (ex1-mid).*([1:N]'/N).^gamma); mid; (mid+ (ex2-mid).*([1:N]'/N).^gamma)];



%% Variance all muscles
Rsquared_slow2= my_Rsquared_coeff(Yslow,hatYslow);
Rsquared_fast2 = my_Rsquared_coeff(Yfast,hatYfast);

Rsquared_slow_unit2 = my_Rsquared_coeff(Yslow,reconstruction_indv(1:168,:));
Rsquared_fast_unit2 = my_Rsquared_coeff(Yfast,reconstruction_indv(169:end,:));

figure()
plot(Rsquared_slow2)
hold on
plot(Rsquared_slow_unit2)
legend('Per leg','Per muscle')

figure
plot(Rsquared_fast2)
hold on
plot(Rsquared_fast_unit2)
legend('Per leg','Per muscle')

% both legs
Rsquared_both= my_Rsquared_coeff([Yslow(:,:);Yfast(:,:)],[hatYslow(:,:);hatYfast(:,:)]);
Rsquared_both_unit = my_Rsquared_coeff([Yslow(:,:);Yfast(:,:)],reconstruction_indv(:,:));

figure
plot(Rsquared_both)
hold on
plot(Rsquared_both_unit)
legend('Per leg','Per muscle')
title('Both leg')
set(gcf,'color','w')
% figure(3)
% ytl={'GLU','TFL','HIP','RF','VL','VM','BF', 'SEMB','SEMT','MG','LG','SOL','PER','TA'};
%%
%% Calculating the R squared per muscle
% mm= 0:12:168;
% mm2=1:12:168;
% for m=1:14
%     Rsquared_slow(m,:) = my_Rsquared_coeff(Yslow(mm2(m):mm(m+1),:),hatYslow(mm2(m):mm(m+1),:));
%     Rsquared_fast(m,:) = my_Rsquared_coeff(Yfast(mm2(m):mm(m+1),:),hatYfast(mm2(m):mm(m+1),:));
%     
%     Rsquared_slow_unit(m,:) = my_Rsquared_coeff(Ymuscles(:,:,m)',Y2asym(:,:,m));
%     Rsquared_fast_unit(m,:) = my_Rsquared_coeff(Ymuscles(:,:,m+14)',Y2asym(:,:,m+14));
%     
% end
% % figure
% muscleOrder={'GLU','TFL','HIP','RF','VL','VM','BF', 'SEMB','SEMT','MG','LG','SOL','PER','TA'};
% ytl= defineMuscleListV2(muscleOrder); %List of muscle
% % ytl={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT', 'VM', 'VL', 'RF', 'HIP','TFL', 'GLU'};
% % ytl(end:-1:1) = ytl(:);
% binw=5;
% 
% mtp=10;
% for m=mtp
%     %     figure
%     %      imagesc(nanmean(Yslow(mm2(m):mm(m+1),481:485),2)',[-1 ,1])
%     %      colormap(flipud(map))
%     figure
%     %     subplot(14,1,m)
%     hold on
%     aux1=conv(Rsquared_slow(m,:),ones(1,binw)/binw,'valid'); %Smoothing
%     plot(aux1,'LineWidth',2,'DisplayName','Per leg','Color',"#0072BD") ;
%     
%     aux1=conv(Rsquared_slow_unit(m,:),ones(1,binw)/binw,'valid'); %Smoothing
%     plot(aux1,'LineWidth',2,'DisplayName','Per muscle','Color',"#A2142F") ;
%     
%     %     plot(movmean(muscles_slow(m,:),5))
%     %     plot(movmean(muscles_slow_unit(m,:),5))
%     %     scatter(1:length(muscles(m,:)),movmean(muscles(m,:),5),'filled')
%     %     legend({'Per leg';'Per muscle'},'AutoUpdate','off');
%     legend('Location','NorthEastOutside','AutoUpdate','off')
%     ylabel(ytl{m+14})
%     %     yline(nanmean(muscles_slow(m,10:30)))
%     
%     pp=patch([40 480 480 40],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
%     
%     %     pp=patch([0 440  440 0],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
%     
%     uistack(pp,'bottom')
% end
% set(gcf,'color','w')
% 
% % figure
% for m=mtp
%     %     subplot(14,1,m)
%     figure
%     hold on
%     aux1=conv(Rsquared_fast(m,:),ones(1,binw)/binw,'valid'); %Smoothing
%     plot(aux1,'LineWidth',2,'DisplayName','Per leg','Color',"#0072BD") ;
%     
%     aux1=conv(Rsquared_fast_unit(m,:),ones(1,binw)/binw,'valid'); %Smoothing
%     plot(aux1,'LineWidth',2,'DisplayName','Per muscle','Color',"#A2142F") ;
%     %         plot(movmean(muscles_fast(m,:),5))
%     %         plot(movmean(muscles_fast_unit(m,:),5))
%     %     scatter(1:length(muscles(m,:)),movmean(muscles(m,:),5),'filled')
%     %     legend({'Per leg';'Per muscle'},'AutoUpdate','off');
%     
%     legend('Location','NorthEastOutside','AutoUpdate','off')
%     ylabel(ytl{m})
%     %     yline(nanmean(muscles_fast(m,10:30)))
%     pp=patch([40 480 480 40],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
%     %     pp=patch([0 440  440 0],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
%     
%     uistack(pp,'bottom')
% end
% set(gcf,'color','w')