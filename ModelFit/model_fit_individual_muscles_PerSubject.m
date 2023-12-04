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

c_constant=0; %1 if you want to use the same regressors for all the participants
plot_visual=1; % 1 if you want to see the results
plotindv=0; % 1 plot individual subjects
indv=1; % 1 individual subjects analsyis; 0 group analysis
adapt=0; % 1 you are looking at data during adaptation; 0 post-adaptation
groupID='NTS' %ID of the group of interest
savedata=1

if contains(groupID,'BAT')
    
    if adapt==1 % For adaptation we are pooling together the data of all 24 participants
        
        fname='dynamicsData_BAT_subj_24_Session_1_RemoveBadMuscles_1_07-August-2023_Adaptation.h5'
    else %Post-adaptation we do it by condiition
        %         fname=['dynamicsData_',groupID,'_subj_12_RemoveBadMuscles_1_PostAdaptation.h5']
        fname=['dynamicsData_',groupID,'_subj_12_Session_1_RemoveBadMuscles_1_01-December-2023_Post-Adaptation.h5']
    end
    
elseif contains(groupID,'CTS') ||  contains(groupID,'NTS')  %TO DO: Add conditions to look at the adaptation
    
    %     fname=['dynamicsData_',groupID,'_subj_5_RemoveBadMuscles_1_PostAdaptation.h5']
    fname=['dynamicsData_',groupID,'_subj_5_Session_1_RemoveBadMuscles_1_04-December-2023_Post-Adaptation.h5']
    
elseif contains(groupID,'C3')
    
    fname= 'dynamicsData_C3_subj_7_Session_1_RemoveBadMuscles_1_17-November-2023_Post-Adaptation.h5'
    
elseif contains(groupID,'AUF')
    if contains(groupID,'V04')
        fname = '21';
    else
        fname = '22';
    end
    if adapt
        fname = ['dynamicsData_AUF_subj_' fname '_RemoveBadMuscles_1_Adapt_TMBase.h5'];
    else
        fname = ['dynamicsData_AUF_subj_' fname '_RemoveBadMuscles_1_OGPost_OGBase.h5']; %h5 file name
    end
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
colorOrder=[ p_orange; p_fade_green; p_fade_blue; p_plum; p_green; p_blue; p_fade_red; p_lime;...
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
        elseif contains(groupID,'CTS') || contains(groupID,'NTS')
            reactive=find(strcmp(epochOfInterest,'SplitPos')==1);
            context= find(strcmp(epochOfInterest,'Adaptation')==1);
            reactive2=find(strcmp(epochOfInterest,'SplitNeg')==1);
            tmbase=find(strcmp(epochOfInterest,'OGNimbus')==1);
        elseif contains(groupID,'C3')
            context= find(strcmp(strtrim(epochOfInterest) ,'Adaptation')==1);
            reactive2=find(strcmp(epochOfInterest,'NegShort_{late}')==1);
            tmbase=find(strcmp(epochOfInterest,'TM base')==1);
        elseif contains(groupID,'AUF')
            load AUFV02_22_IndvLegsC18_ShortPertubations_RemovedBadMuscle_0RemoveBias_0.mat
            
        end
        
    else
        
        if indv==1
            Coefficients=Coefficients_indv(:,:,subj)'; %Cinv means the value for each individual subject
        end
        
        %Pick the conditions that we are going to use for the regression
        epochOfInterest=h5read(fname,'/Epochs')';
        epochOfInterest=deblank(epochOfInterest); %Remove random null values added during the saving of the matrix
        
        if contains(groupID,'BAT')
            %             context= find(strcmp(epochOfInterest,'Optimal')==1);
            context= find(strcmp(epochOfInterest,'Adaptation')==1);
            reactive2=find(strcmp(epochOfInterest,'NegShort_{late}')==1);
            tmbase=find(strcmp(epochOfInterest,'TM base')==1);
            reactive=find(strcmp(epochOfInterest,'Ramp')==1);
        elseif contains(groupID,'CTS') || contains(groupID,'NTS')
            reactive=find(strcmp(epochOfInterest,'SplitPos')==1);
            context= find(strcmp(epochOfInterest,'Adaptation')==1);
            reactive2=find(strcmp(epochOfInterest,'SplitNeg')==1);
            tmbase=find(strcmp(epochOfInterest,'OGNimbus')==1);
        elseif contains(groupID,'C3')
            reactive=find(strcmp(epochOfInterest,'PosShort_{late}')==1);
            context= find(strcmp(strtrim(epochOfInterest) ,'Adaptation')==1);
            reactive2=find(strcmp(epochOfInterest,'NegShort_{late}')==1);
            tmbase=find(strcmp(epochOfInterest,'TM base')==1);
        elseif contains(groupID,'AUF')
            reactive=find(strcmp(epochOfInterest,'PosShort_{la}')==1); %last 10 strides excluding the last 5 from Pos short ramp condition (already reached full slipt)
            context= find(strcmp(epochOfInterest,'Adaptation_{SS}')==1);%last 40 excluding last 5 in the split pos 18 condition (last split, 150 strides splits after 9 short exploration)
            reactive2=find(strcmp(epochOfInterest,'NegShort_{la}')==1); %last 10 excluding last 5 in neg short (full neg split w/o ramp)
            tmbase=find(strcmp(epochOfInterest,'TMBase')==1);
            
        end
        
    end
    
    if adapt==1 %Choosing the regressors for the analysis
        C_regressors=[Coefficients(:,reactive) Coefficients(:,context)]; % EMGreactive and EMGcontext
    else
        C_regressors=[Coefficients(:,reactive2) Coefficients(:,context)]; % EMGreactive and EMGcontext
        
    end
    
    if removebaseline==1
        temp =(EMGdata(1:40,:));
        bias=mean(temp,"omitnan"); %Computing the bias
        
        C_regressors=C_regressors-Coefficients(:,tmbase); %Removing bias from the regressors
        EMGmodel=EMGmodel-bias'; %removing bias from the data
        
        %Removing any NaN for the data. This is needed for the NNMF
        %                 EMGmodel= rmmissing(EMGmodel,2);
        %                 Uf=Uf(:,1:size(EMGmodel,2));
    end
    
    %Organizing the data per muscle
    Coefficients_m=reshape(C_regressors',2,12,size(EMGdata,2)/12); %muscles for regressors 2 number of regressors, 12 number of phase of the gait cycle
    EMGobserved=reshape(EMGmodel(:,1:size(EMGmodel,2))',size(EMGmodel,2),12,size(EMGdata,2)/12); %data
    
    %% Linear regression individual muscles
    
    data=[];
    Regressors_indv=[];
    reconstruction_indv=[];
    %     reconstruction_contextual=[];
    %     reconstruction_reactive=[];
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
        
        %         EMGcontextual(:,:,i)= invEMG_m(:,2,i)* W(2,:,i) ; %yhat_contextual
        %         EMGreactive(:,:,i)=  invEMG_m(:,1,i)* W1,:,i) ; %yhat_contextual
        
        W_transpose(:,:,muscle)=W(:,:,muscle)'; % transposing the dynamics vector
        model{muscle}.C=EMG_coef(:,:,muscle); %saving the regressors. Yhis is necesary for the plotting funciton
        Regressors_indv=[Regressors_indv;EMG_coef(:,:,muscle)];  % Concatenating the regressors for each muscles
        
        %         yContextCurr = EMG_coef(:,2,muscle)* W(2,:,muscle);
        %         yReactiveCurr = EMG_coef(:,1,muscle)* W(1,:,muscle);
        
        %         normContextual(:,subj) = vecnorm(yContextCurr);
        %         normReactive(:,subj) = vecnorm( yReactiveCurr);
        
        %         r = abs(nanmin(yReactiveCurr,[],'all')); %reactive min
        %         c = abs(nanmin(yContextCurr,[],'all')); %contextual min
        
        %         r = abs(nanmin(yReactiveCurr)); %reactive min
        %         c = abs(nanmin(yContextCurr)); %contextual min
        
        reconstruction_indv =[ reconstruction_indv ; EMG_estimated(:,:,muscle)]; % Concatenating the data reconstructed and adding the minimium of the reactive and the contextual
        %         reconstruction_contextual=[reconstruction_contextual; yContextCurr];
        %         reconstruction_reactive=[reconstruction_reactive; yReactiveCurr];
        data_shifted =[ data_shifted ; EMGobserved(:,:,muscle)'];
        data =[ data ; EMGobserved(:,:,muscle)'];  % Concatenating the data
        
        
        % Non-negative matrix factorization (NNMF)
        if indv==0 %tyou can have errors if run in individual participants
            %         YA=[EMGobserved(1:70,:,muscle)' EMGobserved(end-70:end,:,muscle)']; % We only give the algortim data during the transitions
            %         swift=abs(min(EMGobserved(:,:,muscle)',[],'all')); %Finding the most negative value
            %         dataNNMF=YA+swift; %Swifting the data becasue NNMF cant have negative values
            %         data2=EMGobserved(:,:,muscle)';%+swift;
            %         [W2,H] = nnmf(dataNNMF,2);
            %         Casym = W2(:,1:2)';
            %         Cinv=pinv(Casym); %Getting the pseudoinverse of the EMGreactive and EMGcontext
            %         Wasym_NNMF = Cinv'*data2; %x= y/C
            %         NNMF_muscle= ( Casym' * Wasym_NNMF) - swift; %yhat = C
            %         aux1= nanmean(Wasym_NNMF,2).*ones(size(Wasym_NNMF));
            %         NNMF_lower= (Casym' * aux1)- swift; %yhat = C
            %
            %        NNMF=[NNMF ; NNMF_muscle];
            %        NNMF2=[NNMF2 ; NNMF_lower];
            
            % if you want to plot the results of the NNMF - Off for now
            
            %        figure(muscle+20)
            %        subplot 211
            %        plot(W(:,1,muscle))
            %        hold on
            %        plot(aux1(1,:))
            %        plot(Wasym_NNMF(1,:))
            %
            %        subplot 212
            %        plot(W(:,2,muscle))
            %        hold on
            %        plot(aux1(2,:))
            %        plot(Wasym_NNMF(2,:))
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
        aftereffects= mean(EMGobserved(41:45,:,muscle),"omitnan")';
        estimated =mean(EMG_estimated(:,41:45,muscle)',"omitnan")';
        %         PCA_upper=mean(rmmissing(PC(:,41:45))')';
        %         PCA_l=mean(rmmissing(PCA_lower(:,41:45))')';
        %         VIF_F(muscle,:)=vif([model{muscle}.C aftereffects]); %Variance inflation
        %         VIF_F(muscle,:)=vif([model{muscle}.C]); %Variance inflation
        R0 = corrcoef(model{muscle}.C);
        temp=diag(inv(R0))';
        
        VIF_F(muscle,subj) = temp(1);
        
        
        R2(muscle,subj) = my_Rsquared_coeff(aftereffects,estimated,1);
        
        %         R2_upper(muscle,1) = my_Rsquared_coeff(aftereffects,PCA_upper,1);
        %         R2_lower(muscle,1) = my_Rsquared_coeff(aftereffects,PCA_l,1);
        
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
            subplot 311
            hold on
            li{subj}=scatter(muscle,squeeze(mean((W_transpose(41:45,1,muscle)),"omitnan")),25,"filled",'MarkerFaceColor', color(subj,:));
            errorbar(muscle,mdl{subj,muscle}.Coefficients.Estimate(1),mdl{subj,muscle}.Coefficients.SE(1),'k')
            
            subplot 312
            hold on
            li2{subj}=scatter(muscle,squeeze(mean((W_transpose(41:45,2,muscle)),"omitnan")),25,"filled",'MarkerFaceColor', color(subj,:));
            errorbar(muscle,mdl{subj,muscle}.Coefficients.Estimate(2),mdl{subj,muscle}.Coefficients.SE(2),'k')
            
            subplot 313
            hold on
            li3{subj}=scatter(muscle,R2(muscle,subj) ,25,"filled",'MarkerFaceColor', color(subj,:));
            
            
            if plotindv==1
                subplot 311
                ylabel({'Reactive';'AU'})
                yline(0)
                set(gca,'XTick',yt,'XTickLabel',ytl,'FontSize',10)
                legend([li{subj}],subID{subj},'location','Best');
                set(gcf,'color','w')
                
                subplot 312
                ylabel({'Contextual';'AU'})
                yline(0)
                set(gca,'XTick',yt,'XTickLabel',ytl,'FontSize',10)
                %                 lgd=legend([li2{:}],subID{:},'location','Best');
                set(gcf,'color','w')
                
                subplot 313
                ylabel({'R^2'})
                yline(0)
                set(gca,'XTick',yt,'XTickLabel',ytl,'FontSize',10)
                %                 lgd=legend([li3{:}],subID{:},'location','Best');
                set(gcf,'color','w')
                ylim([-1 1])
                yline(.5)
            end
        end
        
        reactive_trace(:,muscle,subj)= W_transpose(:,1,muscle);
        contextual_trace(:,muscle,subj)= W_transpose(:,2,muscle);
        
    end
    
    %% This section is for analysis done for the R01 - ignore for now
    %     data_shifted=nanmean(data_shifted(:,41:45),2);
    %     data_shifted = data_shifted(~isnan(data_shifted));
    %
    %     reconstruction_indv=nanmean(reconstruction_indv(:,41:45),2);
    %     reconstruction_indv=reconstruction_indv(~isnan(reconstruction_indv));
    %
    %     reconstruction_contextual=nanmean(reconstruction_contextual(:,41:45),2);
    %     reconstruction_contextual=reconstruction_contextual(~isnan(reconstruction_contextual));
    %
    %     reconstruction_reactive=nanmean(reconstruction_reactive(:,41:45),2);
    %     reconstruction_reactive=reconstruction_reactive(~isnan(reconstruction_reactive));
    %
    %     r2.total(subj,1)=my_Rsquared_coeff(data_shifted,reconstruction_indv,1);
    %     r2.contextual(subj,1)=my_Rsquared_coeff(data_shifted,reconstruction_contextual,1);
    %     r2.reactive(subj,1)=my_Rsquared_coeff(data_shifted,reconstruction_reactive,1);
    %
    %     r2uncenter.total(subj,1)=my_Rsquared_coeff(data_shifted,reconstruction_indv,0);
    %     r2uncenter.contextual(subj,1)=my_Rsquared_coeff(data_shifted,reconstruction_contextual,0);
    %     r2uncenter.reactive(subj,1)=my_Rsquared_coeff(data_shifted,reconstruction_reactive,0);
    % %
    % index=find(impact>5); % fidn the regressor with high colinearity
    % Uf=Uf(:,1:size(temp2,1));
    % set(gcf,'color','w')
    
    %%
    if plot_visual==1 && subj==n_sub(end) && plotindv==0
        subplot 311
        ylabel({'Reactive';'AU'})
        yline(0)
        set(gca,'XTick',yt,'XTickLabel',ytl,'FontSize',10)
        legend([li{:}],subID,'location','Best')
        set(gcf,'color','w')
        
        subplot 312
        ylabel({'Contextual';'AU'})
        yline(0)
        set(gca,'XTick',yt,'XTickLabel',ytl,'FontSize',10)
        legend([li2{:}],subID,'location','Best')
        set(gcf,'color','w')
    end
    %     if plot_visual==1
    %         model2{1}.C=Regressors_indv; %Save regressors for model
    %         model2{1}.Out= reconstruction_indv; %Save reconstruction of the indiivudal muscles to the model format for plotting
    %         model2{1}.NNMF=NNMF;
    %         model2{1}.NNMF_lower=NNMF2;
    %         analysis=0; %Flag asking if you want to run the analysis
    %         isF=0; %Flag for fast leg
    %
    %
    % %             vizSingleModel_IndvLeg_YoungAdults(model2{1},data,Uf,analysis,[],isF)
    %         % legacy_vizSingleModel_FreeModel_ShortAdaptation_IndvLeg(model2{1},data,Uf,analysis,[],isF)
    %     end
end

%% saving data
if adapt==1 && savedata==1
    
    save([groupID, '_adaptation','_',datestr(now,'dd-mmmm-yyyy'),'.mat'], 'reactive_trace','contextual_trace','R2','fname','subID','mdl','VIF_F')
else
    save([groupID, '_post-adaptation','_',datestr(now,'dd-mmmm-yyyy'),'.mat'], 'reactive_trace','contextual_trace','R2','fname','subID','mdl','VIF_F')
    
end

%% Plot group mean and standart error of the recruitment of the reactive and contextual patterns

reconstruction_contextual=squeeze(mean(rmmissing(contextual_trace(41:45,:,:)),1))'; %contextual
reconstruction_reactive=squeeze(mean(rmmissing(reactive_trace(41:45,:,:)),1))'; % reactive

figure
hold on
for m=1:28 %plotting of all the muscles
    
    subplot 211;hold on
    plot(m,mean(reconstruction_reactive(m,:)),'ok')
    errorbar(m,mean(reconstruction_reactive(m,:)),std(reconstruction_reactive(m,:))/sqrt(1),'k')
    
    subplot 212 ;hold on
    plot(m,mean(reconstruction_contextual(m,:)),'ok')
    errorbar(m,mean(reconstruction_contextual(m,:)),std(reconstruction_contextual(m,:))/sqrt(1),'k')
    
end


%% Plotting time course for individual muscle recruitment

for m=11 %pick the muscle that you want
    
    % Pick the data that you want to plot
    W=W_transpose(40:end,:,m);
    
    figure
    subplot(2,1,1)
    hold on
    scatter(1:length(movmean(W(:,1),binwith)), movmean(W(:,1),binwith),'filled','MarkerFaceColor',"#EDB120") %"#77AC30" )%
    
    title(labels(m).Data)
    legend('Negative Perturbation','AutoUpdate','off')
    pp=patch([0 300 300 0],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
    uistack(pp,'bottom')
    yline(0)
    ylabel({'Reactive';'(A.U)'})
    
    xlabel('strides')
    
    if size(W_transpose(:,:,m),2)>=2
        
        subplot(2,1,1)
        
        hold on
        scatter(1:length(movmean(W(:,2),binwith)), movmean(W(:,2),binwith),'filled','MarkerFaceColor'," #00008B")
        
        legend('Contextual','AutoUpdate','off')
        pp=patch([0 300 300 0],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
        uistack(pp,'bottom')
        yline(0)
        ylabel({'Contextual';'(A.U)'})
        xlabel('strides')
        
    end
    set(gcf,'color','w')
end

%%


