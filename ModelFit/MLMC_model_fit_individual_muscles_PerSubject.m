%Model fit per muscle
%This is a script to run the regression analysis of each muscle using
%pre-determine regressiors see: Mariscal 2023 for references

%Author: Dulce Mariscal, 2023
%%
% Before using this code you need: 
% Process your data, lol 
% Get the .h5 files for this use - preProcessingLinearModel.m


%%
%TODO: Made this script self contained such that we do no need to save data
%form the preProcessingLinearModel.m or made preProcessingLinearModel.m a
%function that can be call in this script


%% Adding path

% upload your path
main='/Users/dulcemariscal/Documents/GitHub/'; %Choose your path
addpath(genpath([main, 'Generalization_Regressions'])) % code for the muscle activity regression
addpath(genpath([main,'labTools'])) % Labtools code - this is the spinal cord of the enitre analysis
% addpath(genpath([main,'splitbelt-EMG-adaptation']))
% addpath(genpath([main,'EMG-LTI-SSM']))
% addpath(genpath([main,'matlab-linsys']))

%% Load real data:
clear all;clc;

%% This is just the saved data - Update accrodingly

c_constant=1; %1 if you want to use the same regressors for all the participants
plot_visual=1; % 1 if you want to see the results
plotindv=0; % 1 plot individual subjects
indv=0; % 1 individual subjects analsyis; 0 group analysis
adapt=0; % 1 you are looking at data during adaptation; 0 post-adaptation
groupID='BATS' %ID of the group of interest
savedata=0;
if contains(groupID,'BAT')
    
    if adapt==1 % For adaptation we are pooling together the data of all 24 participants
        
        fname='dynamicsData_BAT_subj_24_Session_1_RemoveBadMuscles_1_07-August-2023_Adaptation.h5'
    else %Post-adaptation we do it by condiition
        fname=['dynamicsData_',groupID,'_subj_12_RemoveBadMuscles_1_PostAdaptation.h5']
    end
    
elseif contains(groupID,'CTS') %TO DO: Add conditions to look at the adaptation
    
    fname=['dynamicsData_',groupID,'_subj_5_RemoveBadMuscles_1_PostAdaptation.h5']
    
elseif contains(groupID,'C3')
    
    fname= 'dynamicsData_C3_subj_6_Session_1_RemoveBadMuscles_1_30-August-2023_Post-Adaptation.h5'
end

%% Organizaing muscle list.
muscleOrder={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT',...
    'VM', 'VL', 'RF', 'HIP','TFL', 'GLU'};
ytl=([strcat('f',muscleOrder) strcat('s',muscleOrder)]);  %List of muscle
ytl(end:-1:1) = ytl(:);
yt=1:length(ytl);
fs=14;

%% Getting the data
filename=h5read(fname,'/EMGdata');
binwith=10;
[GroupEMGdata,~,U,Ubreaks,~,indvEMGdata,labels,Coefficients,Coefficients_indv]=groupDataToMatrixForm_Update(1:size(filename,3),fname,0);


subID=h5read(fname,'/SubID')';
subID=deblank(subID); %Remove random null values added during the saving of the matrix

Uf=[U;ones(size(U))];
removebaseline=1; %Remove bias flag

if indv==1
    n_sub=size(indvEMGdata,3);
else
    n_sub=1;
end

%% %Defining color for the plots
poster_colors;
colorOrder=[  p_orange; p_red;  p_orange; p_fade_green; p_fade_blue; p_plum; p_green; p_blue; p_fade_red; p_lime;...
    p_yellow; [0 0 0];[0 1 1];p_red; p_orange; p_fade_green; p_fade_blue; p_plum; p_green; p_blue;...
    p_fade_red; p_lime; p_yellow; [0 0 0];[0 1 1]];
color=[ colorOrder; colorOrder;colorOrder];
% color=colormap(jet(size(Yindv,3)));

%% prealoccating variable for speed 

invEMG_m= zeros(12,2,28,'double');
EMG_coef= zeros(12,2,28,'double');
EMG_estimated =zeros(12,size(indvEMGdata,1),28,'double');
W=zeros(2,size(indvEMGdata,1),28,'double');
W_transpose=zeros(size(indvEMGdata,1),2,28,'double');
reactive_trace= zeros(size(indvEMGdata,1),28,n_sub,'double');
contextual_trace= zeros(size(indvEMGdata,1),28,n_sub,'double');
VIF_F= zeros(28,2,'double');
correlation= zeros(28,1,'double');

%% Loop for analysis
for subj=1:n_sub
    
    counter(1:2,subj)=0;
    
    if plotindv==1
        figure
    end
    
    if indv==1 %if individual analysis is 1 then we choose each individual subject
        EMGdata=indvEMGdata(:,:,subj);
    else
        EMGdata=GroupEMGdata;
    end
    
    EMGmodel=EMGdata';  % Transpose the EMG data to match equations
    
    if c_constant==1        
        
        %Pick the conditions that we are going to use for the regression
        epochOfInterest=h5read(fname,'/Epochs')';
        epochOfInterest=deblank(epochOfInterest);
        
        if contains(groupID,'BAT')
            load BAT_24_IndvLegsC16_ShortPertubations_RemovedBadMuscle_1.mat
            Coefficients=C;
            context= find(strcmp(epochOfInterest,'Optimal')==1);
            %             context= find(strcmp(epochOfInterest,'Adaptation')==1);
            reactive2=find(strcmp(epochOfInterest,'NegShort_{late}')==1);
            tmbase=find(strcmp(epochOfInterest,'TM base')==1);
            reactive=find(strcmp(epochOfInterest,'Ramp')==1);
        elseif contains(groupID,'CTS')
            reactive=find(strcmp(epochOfInterest,'SplitPos')==1);
            context= find(strcmp(epochOfInterest,'Adaptation')==1);
            reactive2=find(strcmp(epochOfInterest,'SplitNeg')==1);
            tmbase=find(strcmp(epochOfInterest,'OGNimbus')==1);
        elseif contains(groupID,'C3')
            context= find(strcmp(strtrim(epochOfInterest) ,'Adaptation')==1);
            reactive2=find(strcmp(epochOfInterest,'NegShort_{early}')==1);
            tmbase=find(strcmp(epochOfInterest,'TM base')==1);
        end
        
    else
        
        if indv==1
            Coefficients=Coefficients_indv(:,:,subj)'; %Cinv means the value for each individual subject
        end
        
        %Pick the conditions that we are going to use for the regression
        epochOfInterest=h5read(fname,'/Epochs')';
        epochOfInterest=deblank(epochOfInterest); %Remove random null values added during the saving of the matrix
        
        if contains(groupID,'BAT')
            context= find(strcmp(epochOfInterest,'Optimal')==1);
            reactive2=find(strcmp(epochOfInterest,'NegShort_{late}')==1);
            tmbase=find(strcmp(epochOfInterest,'TM base')==1);
            reactive=find(strcmp(epochOfInterest,'Ramp')==1);
        elseif contains(groupID,'CTS')
            reactive=find(strcmp(epochOfInterest,'SplitPos')==1);
            context= find(strcmp(epochOfInterest,'Adaptation')==1);
            reactive2=find(strcmp(epochOfInterest,'SplitNeg')==1);
            tmbase=find(strcmp(epochOfInterest,'OGNimbus')==1);
        elseif contains(groupID,'C3')
            reactive=find(strcmp(epochOfInterest,'PosShort_{late}')==1);
            context= find(strcmp(strtrim(epochOfInterest) ,'Adaptation')==1);
            reactive2=find(strcmp(epochOfInterest,'NegShort_{late}')==1);
            tmbase=find(strcmp(epochOfInterest,'TM base')==1);
            
        end
        
    end
    
    if adapt==1 %Choosing the regressors for the analysis
        C_regressors=[Coefficients(:,reactive) Coefficients(:,context)]; % EMGreactive and EMGcontext
    else
        C_regressors=[Coefficients(:,reactive2) Coefficients(:,context)]; % EMGreactive and EMGcontext
        
    end
    
    if removebaseline==1
        temp = rmmissing(EMGdata(1:40,:));
        bias=mean(temp); %Computing the bias
        
        C_regressors=C_regressors-Coefficients(:,tmbase); %Removing bias from the regressors
        EMGmodel=EMGmodel-bias'; %removing bias from the data
        
        %Removing any NaN for the data. This is needed for the NNMF
%                 EMGmodel= rmmissing(EMGmodel,2);
%                 Uf=Uf(:,1:size(EMGmodel,2));
% if contains(groupID,'BATR')
%          C_regressors=[Coefficients(:,reactive2) Coefficients(:,context) Coefficients(:,tmbase)];
% else
%         C_regressors=[Coefficients(:,reactive2) Coefficients(:,context) bias'];
% end
    end
    
%     %Organizing the data per muscle
    Coefficients_m=reshape(C_regressors',2,12,size(EMGdata,2)/12); %muscles for regressors 2 number of regressors, 12 number of phase of the gait cycle
    EMGobserved=reshape(EMGmodel(:,1:size(EMGmodel,2))',size(EMGmodel,2),12,size(EMGdata,2)/12); %data
    
    %% Linear regression individual muscles
    
    data=[];
    Regressors_indv=[];
    reconstruction_indv=[];
    data_shifted =[];
    NNMF=[]; NNMF2=[];
    
    for muscle=1:28
%         yContextCurr=[];
%         yReactiveCurr=[];
        
        EMG_coef(:,:,muscle)=Coefficients_m(:,:,muscle)'./vecnorm(Coefficients_m(:,:,muscle)'); %Getting the unit vector of the regressors
        
        if all(isnan(EMG_coef(:,:,muscle)))
            invEMG_m(:,:,muscle)=(EMG_coef(:,:,muscle)); %if nan can't get inverse, carry the nan alone.
        else
            invEMG_m(:,:,muscle)=pinv(EMG_coef(:,:,muscle)'); % Getting the inverse, if data is not nan.
        end
        
        W(:,:,muscle) =invEMG_m(:,:,muscle)'*EMGobserved(:,:,muscle)'; % W = EMG^-1 * EMGdata^m
        EMG_estimated(:,:,muscle)=  EMG_coef(:,:,muscle)* W(:,:,muscle) ; %EMGhat = estimated muscle activity
        

        
        W_transpose(:,:,muscle)=W(:,:,muscle)'; % transposing the dynamics vector
        model{muscle}.C=EMG_coef(:,:,muscle); %saving the regressors. Yhis is necesary for the plotting funciton
        Regressors_indv=[Regressors_indv;EMG_coef(:,:,muscle)];  % Concatenating the regressors for each muscles
        

        
        reconstruction_indv =[ reconstruction_indv ;EMG_estimated(:,:,muscle)]; % Concatenating the data reconstructed and adding the minimium of the reactive and the contextual
        data_shifted =[ data_shifted ; EMGobserved(:,:,muscle)'];
        data =[ data ; EMGobserved(:,:,muscle)'];  % Concatenating the data
        
        %% PCA
        if indv==0 %tyou can have errors if run in individual participants
            YA=[EMGobserved(:,:,muscle)']; % We only give the algortim data during the transitions
            
            [pp,cc_,aa_]=pca(YA');
            Wpca=[pp(:,1:2)];
            Winv=pinv(Wpca);
            dynamics= EMGobserved(:,:,muscle)*Winv';
            PC= Wpca * dynamics' ;
            aux1= nanmean(dynamics).*ones(size(dynamics));
            PCA_lower= (Wpca * aux1'); %yhat = C
      
        end
        
       %%
        % Checking for colinearity and correlation between the regresssions
        R=corrcoef(model{muscle}.C); % correlation coeficient
        correlation(muscle,1)=R(2);
        aftereffects= mean(rmmissing(EMGobserved(41:45,:,muscle)))';
        estimated =mean(rmmissing(EMG_estimated(:,41:45,muscle))')';
        PCA_upper=mean(rmmissing(PC(:,41:45))')';
        PCA_l=mean(rmmissing(PCA_lower(:,41:45))')';
%         VIF_F(muscle,:)=vif([model{muscle}.C aftereffects]); %Variance inflation
%         VIF_F(muscle,:)=vif([model{muscle}.C]); %Variance inflation
        R0 = corrcoef(model{muscle}.C);
        VIF_F(muscle,:) = diag(inv(R0))';
        
        
        R2(muscle,1) = my_Rsquared_coeff(aftereffects,estimated,1);
        R2_upper(muscle,1) = my_Rsquared_coeff(aftereffects,PCA_upper,1);
        R2_lower(muscle,1) = my_Rsquared_coeff(aftereffects,PCA_l,1);
        
        %fit linear model for post-adaptation - the intercept is off. We
        %can use this to get the CI of the data. 
       
        %%TODO: Check why we have the same values for both regressors
       
        mdl{subj,muscle}= fitlm(EMG_coef(:,:,muscle),aftereffects,'VarNames',{'Reactive','Contextual', labels(muscle).Data(1:end-1)},'Intercept',false);
        
        if mdl{subj,muscle}.Coefficients.Estimate(2)+mdl{subj,muscle}.Coefficients.SE(2)>0 &&  mdl{subj,muscle}.Coefficients.Estimate(2)-mdl{subj,muscle}.Coefficients.SE(2)<0
            counter(2,subj)=counter(2,subj)+1;
        end
        
        if mdl{subj,muscle}.Coefficients.Estimate(1)+mdl{subj,muscle}.Coefficients.SE(1)>0 &&  mdl{subj,muscle}.Coefficients.Estimate(1)-mdl{subj,muscle}.Coefficients.SE(1)<0
            counter(1,subj)=counter(1,subj)+1;
        end
        
        if plot_visual==1
%             figure(1)
%             subplot 211
%             hold on
%             li{subj}=scatter(muscle,squeeze(mean(rmmissing(W_transpose(41:45,1,muscle)))),25,"filled",'MarkerFaceColor', color(subj,:));
%             errorbar(muscle,mdl{subj,muscle}.Coefficients.Estimate(1),mdl{subj,muscle}.Coefficients.SE(1),'k')
%             
%             subplot 212
%             hold on
%             li2{subj}=scatter(muscle,squeeze(mean(rmmissing(W_transpose(41:45,2,muscle)))),25,"filled",'MarkerFaceColor', color(subj,:));
%             errorbar(muscle,mdl{subj,muscle}.Coefficients.Estimate(2),mdl{subj,muscle}.Coefficients.SE(2),'k')
%             
            
            if plotindv==1
                subplot 211
                ylabel({'Reactive';'AU'})
                yline(0)
                set(gca,'XTick',yt,'XTickLabel',ytl,'FontSize',10)
                legend([li{subj}],subID{subj},'location','Best')
                set(gcf,'color','w')
                
                subplot 212
                ylabel({'Contextual';'AU'})
                yline(0)
                set(gca,'XTick',yt,'XTickLabel',ytl,'FontSize',10)
                legend([li2{subj}],subID{subj},'location','Best')
                set(gcf,'color','w')
            end
            
            figure(3)
            subplot 211
            hold on
            li{subj}=scatter(muscle,R2(muscle,1),100,"filled",'MarkerFaceColor',  [0 0 0]);
%             scatter(muscle,R2_upper(muscle,1),"_",'MarkerEdgeColor', 'r','LineWidth',5);
            
            subplot 212
            hold on
            li2{subj}=scatter(muscle,R2(muscle,1),100,"filled",'MarkerFaceColor',  [0 0 0]);
%             scatter(muscle,R2_upper(muscle,1),"_",'MarkerEdgeColor', 'r','LineWidth',5);
           
        end
        
        reactive_trace(:,muscle,subj)= W_transpose(:,1,muscle);
        contextual_trace(:,muscle,subj)= W_transpose(:,2,muscle);
        
        if muscle==100
            
            model2{2}.C=EMG_coef(:,:,muscle); %Save regressors for model
            model2{2}.Out= EMG_estimated(:,:,muscle); %Save reconstruction of the indiivudal muscles to the model format for plotting
            data2= EMGobserved(:,:,muscle)';
            analysis=0; %Flag asking if you want to run the analysis
            isF=0; %Flag for fast leg
            
            
            vizSingleModel_IndvLeg_YoungAdults(model2{2},data2,Uf,analysis,labels(muscle).Data,isF)
            
        end
        
    end
    
    
    %%
    if plot_visual==1 && subj==n_sub(end) && plotindv==0
        subplot 211
        ylabel({'Reactive';'AU'})
        yline(0)
        set(gca,'XTick',yt,'XTickLabel',ytl,'FontSize',10)
        legend([li{:}],subID,'location','Best')
        set(gcf,'color','w')
        
        subplot 212
        ylabel({'Contextual';'AU'})
        yline(0)
        set(gca,'XTick',yt,'XTickLabel',ytl,'FontSize',10)
        legend([li2{:}],subID,'location','Best')
        set(gcf,'color','w')
    end
    if plot_visual==1
        model2{1}.C=Regressors_indv; %Save regressors for model
        model2{1}.Out= reconstruction_indv; %Save reconstruction of the indiivudal muscles to the model format for plotting
%         model2{1}.NNMF=NNMF;
%         model2{1}.NNMF_lower=NNMF2;
        analysis=0; %Flag asking if you want to run the analysis
        isF=0; %Flag for fast leg
        
        
            vizSingleModel_IndvLeg_YoungAdults(model2{1},data,Uf,analysis,[],isF)


    end
end

%% saving data
if adapt==1 && savedata==1
    
    save([groupID, '_adaptation.mat'], 'W','W_transpose')
else
    save([groupID, '_post-adaptation.mat'], 'W','W_transpose')
    
end

%% Plotting time course for individual muscle recruitment 

for m=10%pick the muscle that you want
    
    % Pick the data that you want to plot
    W=W_transpose(40:end,:,m);
   
    figure
    subplot(2,1,1)
    hold on
    scatter(1:length(movmean(W(:,1),binwith)), -movmean(W(:,1),binwith),'filled','MarkerFaceColor',"#EDB120") %"#77AC30" )%
    
    title(labels(m).Data)
    legend('Negative Perturbation','AutoUpdate','off')
%     pp=patch([0 300 300 0],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
%     uistack(pp,'bottom')
%     yline(0)
%     ylabel({'Reactive';'(A.U)'})

    xlabel('strides')
    
    if size(W_transpose(:,:,m),2)>=2
  
        subplot(2,1,1)
        
        hold on
        scatter(1:length(movmean(W(:,2),binwith)), movmean(W(:,2),binwith),'filled','MarkerFaceColor'," #00008B")
        
%                 hold on
%         scatter(1:length(movmean(W(:,3),binwith)), movmean(W(:,3),binwith),'filled','MarkerFaceColor',"r")
        
        legend('Contextual','AutoUpdate','off')  
%         pp=patch([0 300 300 0],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
%         uistack(pp,'bottom')
        yline(0)
        ylabel({'W^{fSEMT}';'(A.U)'})
        xlabel('strides')
    
    end
    set(gcf,'color','w') 
end

%% Plotting the scatter plot for post-adaptation 
% I saved the W from the BATS and BATR groups post-adaptation 

load('musclesLabels.mat')

load('BATS_post-adaptation.mat')
OG=W;

load('BATR_post-adaptation.mat')
TM=W;


steps=41:45;

hold on
for c=1:2
    figure
    hold on
    if c==1
        TM(1,:,:)=TM(1,:,:);
        OG(1,:,:)=OG(1,:,:);
    end
    for i=1:28
        
        if i<=3 %'sHIP'
            Li{1}=scatter(nanmean(TM(c,steps,i)),nanmean(OG(c,steps,i)),100,"filled",'MarkerFaceColor',"white","MarkerEdgeColor", "#77AC30");
            text(nanmean(TM(c,steps,i)),nanmean(OG(c,steps,i)),{labels(i).Data(1:end-1)})
            
        elseif i>=4 &&  i<=9 %sTHIGH
            Li{2}=scatter(nanmean(TM(c,steps,i)),nanmean(OG(c,steps,i)),100,"filled",'MarkerFaceColor',"white","MarkerEdgeColor","#7E2F8E");
            text(nanmean(TM(c,steps,i)),nanmean(OG(c,steps,i)),{labels(i).Data(1:end-1)})
            
        elseif i>=10 &&  i<=14 %sSHANK
            Li{3}=scatter(nanmean(TM(c,steps,i)),nanmean(OG(c,steps,i)),100,"filled",'MarkerFaceColor',"white","MarkerEdgeColor","#D95319");
            text(nanmean(TM(c,steps,i)),nanmean(OG(c,steps,i)),{labels(i).Data(1:end-1)})
            
        elseif i>=15 &&  i<=17 %fHIP
            Li{4}=scatter(nanmean(TM(c,steps,i)),nanmean(OG(c,steps,i)),100,"filled",'MarkerFaceColor', "#77AC30");
            text(nanmean(TM(c,steps,i)),nanmean(OG(c,steps,i)),{labels(i).Data(1:end-1)})
            
        elseif i>=18 &&  i<=23 %fTHIGH
            Li{5}=scatter(nanmean(TM(c,steps,i)),nanmean(OG(c,steps,i)),100,"filled",'MarkerFaceColor', "#7E2F8E");
            text(nanmean(TM(c,steps,i)),nanmean(OG(c,steps,i)),{labels(i).Data(1:end-1)})
            
        elseif i>=24 %fSHANK
            Li{6}=scatter(nanmean(TM(c,steps,i)),nanmean(OG(c,steps,i)),100,"filled",'MarkerFaceColor', "#D95319");
            text(nanmean(TM(c,steps,i)),nanmean(OG(c,steps,i)),{labels(i).Data(1:end-1)})
        end
        
    end
    
    if c==1
        xx=-.5:0.1:1.5;
        plot(xx,xx,'k')
        
    ylabel({'OG_{Reactive}';'A.U'})
    xlabel({'TM_{Reactive}';'A.U'})
    elseif c==2
        yline(0)
        xline(0)
        xlim([-0.5 1])
        ylim([-0.5 1])
        
    ylabel({'OG_{contextual}';'A.U'})
    xlabel({'TM_{contextual}';'A.U'})
    end
    set(gcf,'color','w')
    set(gcf,'renderer','painters')
    axis square
    legend([Li{:}],{'sHIP';'sTHIGH';'sSHANK';'fHIP';'fTHIGH';'fSHANK'})
end





%% Time courses for the post-adaptation data. 
% I saved the W from the BATS and BATR groups post-adaptation 

load('musclesLabels.mat')
load('BATR_post-adaptation.mat')
TM=W;

load('BATS_post-adaptation.mat')
OG=W;
steps=41:45;

 binwith=5;
 analysis=0;
 trace=2;
 sz = 50;
 
 facecolor=["#EDB120";" #00008B"];
 figure
 for c=1:2
     
     if c==1
         TM(1,:,:)=-TM(1,:,:);
         OG(1,:,:)=-OG(1,:,:);
     end
     
     for i=25
         
         
         data=[TM(c,41:end,i)];
         
         subplot(2,1,c)
         hold on
         
         L{1}=scatter(1:length(movmean(data,binwith)), movmean(data,binwith),sz,'filled','MarkerFaceColor', facecolor(1));

         if size(TM(:,:,i),2)>=2
             
             subplot(2,1,c)
             data=[OG(c,41:end,i)];
             hold on
             
             L{2} = scatter(1:length(movmean(data,binwith)), movmean(data,binwith),sz,'filled','MarkerFaceColor',"white","MarkerEdgeColor",facecolor(2));

         end
         
         title(labels(i).Data)
         set(gcf,'color','w')
         legend([L{:}],{'Treadmill','Overground'},'AutoUpdate','off')
         yline(0)
         
         xlabel('strides')
         
         if c==1
             ylabel({'Reactive';'(A.U)'})
         else
             ylabel({'Contextual';'(A.U)'})
         end
         
     end
     
     xlim([0 200])
 end