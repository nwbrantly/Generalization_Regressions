%Model fit per muscle

%% Adding path

% upload your path
main='/Users/dulcemariscal/Documents/GitHub/';
addpath(genpath([main, 'Generalization_Regressions']))
addpath(genpath([main,'labTools']))
addpath(genpath([main,'splitbelt-EMG-adaptation']))
addpath(genpath([main,'EMG-LTI-SSM']))
addpath(genpath([main,'matlab-linsys']))

%% Load real data:
clear all;clc
%% Free model - Linear regression using psedoinverse
%% This is just the saved data - Update accrodingly
% fname = ['dynamicsData_', groupID, '_subj_3_RemoveBadMuscles_0.h5']

c_constant=1;
plot_visual=1;
plotindv=0;
indv=0;
adapt=0;
groupID='BATR'

if contains(groupID,'BAT')

    if adapt==1
        fname=['dynamicsData_BAT_subj_24_Session_1_RemoveBadMuscles_1_07-August-2023_Adaptation.h5']
    else
        fname=['dynamicsData_',groupID,'_subj_12_RemoveBadMuscles_1_PostAdaptation.h5']
    end
       
elseif contains(groupID,'CTS')
   
    fname=['dynamicsData_',groupID,'_subj_5_RemoveBadMuscles_1_PostAdaptation.h5']

elseif contains(groupID,'C3')
    
    fname=['dynamicsData_C3_subj_6_Session_1_RemoveBadMuscles_1_30-August-2023_Post-Adaptation.h5']
end
    
%% 
muscleOrder={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT',...
    'VM', 'VL', 'RF', 'HIP','TFL', 'GLU'};
ytl=([strcat('f',muscleOrder) strcat('s',muscleOrder)]);  %List of muscle
ytl(end:-1:1) = ytl(:);
yt=1:length(ytl);
fs=14;

%% Getting the data
EMGdata=h5read(fname,'/EMGdata');
binwith=10;
[Y,~,U,Ubreaks,Ysum,Yindv,labels,C,Cinv]=groupDataToMatrixForm_Update(1:size(EMGdata,3),fname,0);

% color=colormap(jet(size(Yindv,3)));
poster_colors;
colorOrder=[p_red; p_orange; p_fade_green; p_fade_blue; p_plum; p_green; p_blue; p_fade_red; p_lime;...
    p_yellow; [0 0 0];[0 1 1];p_red; p_orange; p_fade_green; p_fade_blue; p_plum; p_green; p_blue;...
    p_fade_red; p_lime; p_yellow; [0 0 0];[0 1 1]];
color=[ colorOrder; colorOrder;colorOrder];

%% Organizing the data

subID=h5read(fname,'/SubID')';
subID=deblank(subID); %Remove random null values added during the saving of the matrix

Uf=[U;ones(size(U))];
removebaseline=1; %Remove bias flag

if indv==1
    n_sub=size(Yindv,3);
else
     n_sub=1;
end

for subj=1:n_sub
    
   counter(1:2,subj)=0;
    
    if plotindv==1
        figure
    end
            
    if indv==1
        Yasym=Yindv(:,:,subj);
    else
        Yasym=Y;
    end 
    
    Ymodel=Yasym';  % Transpose the EMG data to match equations
    
    if c_constant==1
        %Pick the conditions that we are going to use for the regression
        load BAT_24_IndvLegsC16_ShortPertubations_RemovedBadMuscle_1.mat
%         C=C';



        %Pick the conditions that we are going to use for the regression
%         epochOfInterest=h5read(fname,'/Epochs')';
%         epochOfInterest=deblank(epochOfInterest);
        
        if contains(groupID,'BAT')
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
            C=Cinv(:,:,subj)';
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
    
    if adapt==1
        
        Casym=[C(:,reactive) C(:,context)]; % EMGreactive and EMGcontext
    else
        Casym=[C(:,reactive2) C(:,context)]; % EMGreactive and EMGcontext
        
    end
    
    if removebaseline==1
        bias=nanmean(Yasym(1:40,:)); %Computing the bias
        
        Casym=Casym-C(:,tmbase); %Removing bias from the regressors
        Ymodel=Ymodel-bias'; %removing bias from the data
        
        %Removing any NaN for the data
%         Ymodel= rmmissing(Ymodel,2);
%         Uf=Uf(:,1:size(Ymodel,2));
    end
    
    %Organizing the data per muscle
    Cmuscles=reshape(Casym',2,12,size(Yasym,2)/12); %muscles for regressors 2 number of regressors, 12 number of phase of the gait cycle
    Ymuscles=reshape(Ymodel(:,1:size(Ymodel,2))',size(Ymodel,2),12,size(Yasym,2)/12); %data
    
    %% Linear regression individual muscles
    
    data=[];
    C_indv=[];
    reconstruction_indv=[];
    reconstruction_contextual=[];
    reconstruction_reactive=[];
    data_shifted =[];
    NNMF=[]; NNMF2=[];
        
    for muscles=1:28
        yContextCurr=[];
        yReactiveCurr=[];
        
        unit(:,:,muscles)=Cmuscles(:,:,muscles)'./vecnorm(Cmuscles(:,:,muscles)'); %Getting the unit vector of the regressors
        %         temp(:,:,i)=pinv(unit(:,:,i)'); % Getting the inverse

        if all(isnan(unit(:,:,muscles)))
            temp(:,:,muscles)=(unit(:,:,muscles)); %if nan can't get inverse, carry the nan alone.
        else
            temp(:,:,muscles)=pinv(unit(:,:,muscles)'); % Getting the inverse, if data is not nan.
        end
        
        X2asym(:,:,muscles) =temp(:,:,muscles)'*Ymuscles(:,:,muscles)'; %x= y/C
        Y2asym(:,:,muscles)=  unit(:,:,muscles)* X2asym(:,:,muscles) ; %yhat
        
%         Ycontextual(:,:,i)=  unit(:,2,i)* X2asym(2,:,i) ; %yhat_contextual
%         Yreactive(:,:,i)=  unit(:,1,i)* X2asym(1,:,i) ; %yhat_contextual
        
        temp2(:,:,muscles)=X2asym(:,:,muscles)'; % transposing the dynamics vector
        model{muscles}.C=unit(:,:,muscles); %saving the regressors. Yhis is necesary for the plotting funciton
        C_indv=[C_indv;unit(:,:,muscles)];  % Concatenating the regressors for each muscles
        
        yContextCurr = unit(:,2,muscles)* X2asym(2,:,muscles);
        yReactiveCurr = unit(:,1,muscles)* X2asym(1,:,muscles);
        
        normContextual(:,subj) = vecnorm(yContextCurr);
        normReactive(:,subj) = vecnorm( yReactiveCurr);
             
%         r = abs(nanmin(yReactiveCurr,[],'all')); %reactive min
%         c = abs(nanmin(yContextCurr,[],'all')); %contextual min

     r = abs(nanmin(yReactiveCurr)); %reactive min
    c = abs(nanmin(yContextCurr)); %contextual min       
        
        reconstruction_indv =[ reconstruction_indv ; Y2asym(:,:,muscles)]; % Concatenating the data reconstructed and adding the minimium of the reactive and the contextual
        reconstruction_contextual=[reconstruction_contextual; yContextCurr];
        reconstruction_reactive=[reconstruction_reactive; yReactiveCurr];
        data_shifted =[ data_shifted ; Ymuscles(:,:,muscles)'];
        data =[ data ; Ymuscles(:,:,muscles)'];  % Concatenating the data


        % NNMF
%         YA=[Ymuscles(1:70,:,muscles)' Ymuscles(end-70:end,:,muscles)'];
% %         YA=[Ymuscles(:,:,muscles)''];
%         swift=abs(min(Ymuscles(:,:,muscles)',[],'all'));
%         dataNNMF=YA+swift;
%         data2=Ymuscles(:,:,muscles)';%+swift;
%         [W,H] = nnmf(dataNNMF',2);
% %         [W,H] = pca(dataNNMF);
%         Casym = W(:,1:2)';
%         Cinv=pinv(Casym); %Gettign the pseudoinverse of the EMGreactive and EMGcontext
%         Wasym_NNMF = Cinv'*data2; %x= y/C
%         NNMF_muscle= ( Casym' * Wasym_NNMF) - swift; %yhat = C
%         aux1= nanmean(Wasym_NNMF,2).*ones(size(Wasym_NNMF));
%         NNMF_lower= (Casym' * aux1)- swift; %yhat = C
%         
%        NNMF=[NNMF ; NNMF_muscle];
%        NNMF2=[NNMF2 ; NNMF_lower];

%        figure(muscles+20)
%        subplot 211
%        plot(temp2(:,1,muscles))
%        hold on
%        plot(aux1(1,:))
%        plot(Wasym_NNMF(1,:)) 
%        
%        subplot 212
%        plot(temp2(:,2,muscles))
%        hold on
%        plot(aux1(2,:))
%         plot(Wasym_NNMF(2,:))
        
        % Checking for colinearity and correlation between the regresssions
        temp5=corrcoef(model{muscles}.C); % correlation coeficient
        correlation(muscles,1)=temp5(2);
        temp10(muscles,:)=vif([model{muscles}.C nanmean(Ymuscles(41:51,:,muscles))']); %Variance inflation
        temp5=vif([model{muscles}.C]);
        impact(muscles,:)=temp5;
        
        %fit linear model
        mdl{subj,muscles}= fitlm(unit(:,:,muscles),nanmean(Ymuscles(41:45,:,muscles))','VarNames',{'Reactive','Contextual', labels(muscles).Data(1:end-1)},'Intercept',false);
        mdl2{subj,muscles}= fitlm(unit(:,:,muscles),nanmean(Ymuscles(41:45,:,muscles))','VarNames',{'Reactive','Contextual', labels(muscles).Data(1:end-1)},'Intercept',true);
        
        if mdl{subj,muscles}.Coefficients.Estimate(2)+mdl{subj,muscles}.Coefficients.SE(2)>0 &&  mdl{subj,muscles}.Coefficients.Estimate(2)-mdl{subj,muscles}.Coefficients.SE(2)<0
            counter(2,subj)=counter(2,subj)+1;
        end
        
        if mdl{subj,muscles}.Coefficients.Estimate(1)+mdl{subj,muscles}.Coefficients.SE(1)>0 &&  mdl{subj,muscles}.Coefficients.Estimate(1)-mdl{subj,muscles}.Coefficients.SE(1)<0
            counter(1,subj)=counter(1,subj)+1;
        end
        
        if plot_visual==1 
            figure(1)
            subplot 211
            hold on
            li{subj}=scatter(muscles,squeeze(nanmean(temp2(41:45,1,muscles))),25,"filled",'MarkerFaceColor', color(2,:));
            errorbar(muscles,mdl{subj,muscles}.Coefficients.Estimate(1),mdl{subj,muscles}.Coefficients.SE(1),'k')
            
            subplot 212
            hold on
            li2{subj}=scatter(muscles,squeeze(nanmean(temp2(41:45,2,muscles))),25,"filled",'MarkerFaceColor', color(2,:));
            errorbar(muscles,mdl{subj,muscles}.Coefficients.Estimate(2),mdl{subj,muscles}.Coefficients.SE(2),'k')
            
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
        end
        
        reactive_trace(:,muscles,subj)= temp2(:,1,muscles);
        contextual_trace(:,muscles,subj)= temp2(:,2,muscles);

    end
    
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
    
    if plot_visual==1 && subj==n_sub(end) && plotindv==0
%         figure
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
    model2{1}.C=C_indv; %Save regressors for model
    model2{1}.Out= reconstruction_indv; %Save reconstruction of the indiivudal muscles to the model format for plotting
    model2{1}.NNMF=NNMF;
    model2{1}.NNMF_lower=NNMF2;
    analysis=0; %Flag asking if you want to run the analysis
    isF=0; %Flag for fast leg
    

%     vizSingleModel_IndvLeg_YoungAdults(model2{1},data,Uf,analysis,[],isF)
    % legacy_vizSingleModel_FreeModel_ShortAdaptation_IndvLeg(model2{1},data,Uf,analysis,[],isF)
  end
end


%%
    reconstruction_contextual=squeeze(nanmean(reactive_trace(41:45,:,:),1))' ;  
    reconstruction_reactive=squeeze(nanmean(reactive_trace(41:45,:,:),1))'; 
%%
    figure 
    hold on
    for m=1:28
        subplot 211;hold on
        plot(m,nanmean(reconstruction_reactive(m,:)),'ok')
        errorbar(m,nanmean(reconstruction_reactive(m,:)),nanstd(reconstruction_reactive(m,:))/sqrt(6),'k')
        
        subplot 212 ;hold on
        plot(m,nanmean(reconstruction_contextual(m,:)),'ok')
        errorbar(m,nanmean(reconstruction_contextual(m,:)),nanstd(reconstruction_contextual(m,:))/sqrt(6),'k')
    
    end
    
    subplot 211;hold on
    yline(0)
    ylabel('W_{reactive}')
    set(gca,'XTick',yt,'XTickLabel',ytl,'FontSize',10)
    subplot 212;hold on
    yline(0)
    ylabel('W_{contextual}')
    set(gca,'XTick',yt,'XTickLabel',ytl,'FontSize',10)
%%
% SLA=h5read(fname,'/SLA')';
% SLA_post1= nanmean(SLA(:, 41:45),2);
% 
% figure
% scatter(counter(2,:),SLA_post1,'filled')
% disp('# contextual vs SLA')
% [rho,pval] = corr(counter(2,:)',SLA_post1)
% xlabel('# contextual ')
% ylabel('SLA')
% 
% figure
% scatter(counter(1,:),counter(2,:),'filled')
% xlabel('# reactive')
% ylabel('# contextual')
% disp('# reactive vs # contextual')
% [rho,pval] = corr(counter(1,:)',counter(2,:)')
% 
% %%
% 
% After_norm_contextual= nanmean(normContextual(41:45,:),1);
% 
% After_norm_reactive= nanmean(normReactive(41:45,:),1);
% 
% figure 
% scatter(After_norm_contextual,After_norm_reactive,'filled')
% ylabel('# reactive')
% xlabel('# contextual')
% disp('Norm contextual vs Norm Reactive')
% [rho,pval] = corr(After_norm_contextual',After_norm_reactive')
% 
% 
% %%
% 
% 
% figure 
% scatter(After_norm_contextual,SLA_post1,'filled')
% ylabel('SLA')
% xlabel('Norm of contextual')
% disp('Norm contextual vs SLA')
% [rho,pval] = corr(After_norm_contextual',SLA_post1)
% 
% 
% 
% figure 
% scatter(After_norm_reactive,SLA_post1,'filled')
% ylabel('SLA')
% xlabel('Norm of reactive')
% disp('Norm contextual vs SLA')
% [rho,pval] = corr(After_norm_reactive',SLA_post1)
% 
% 
% %%
% % save([groupID,'_Indv_',num2str(indv),'_Adapt_',num2str(adapt),'_regressorConstant_',num2str(c_constant),'_removeBias_',num2str(removebaseline),'.mat'],'subID','ytl','reactive_trace','contextual_trace','fname','r2','r2uncenter','mdl')
% 
% 
% 
% %%
% % Uf=Uf(:,1:size(temp2,1));
% % %%
% % 
% % 
% % 
% % %%  Plotting
% % model2{1}.C=C_indv; %Save regressors for model
% % model2{1}.Out= reconstruction_indv; %Save reconstruction of the indiivudal muscles to the model format for plotting
% % analysis=0; %Flag asking if you want to run the analysis
% % isF=0; %Flag for fast leg
% % 
% % 
% % legacy_vizSingleModel_FreeModel_ShortAdaptation_IndvLeg2(model2{1},data,Uf,analysis,[],isF)
% % % legacy_vizSingleModel_FreeModel_ShortAdaptation_IndvLeg(model2{1},data,Uf,analysis,[],isF)
% % % save('BATS_indv_muscles.mat','X2asym')
% % % end
% % 
%% Plotting time course for individual muscles
analysis=0

for i=25
    %pick the muscle that you want
    
    % Pick the data that you want to plot
    %    Xasym=[temp2(1:40,:,i);nan(1,size(temp2,2));temp2(41:480,:,i);nan(1,size(temp2,2));temp2(481:end,:,i)];
    Xasym=[temp2(40:end,:,i)];
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
%     legacy_vizSingleModel_FreeModel_ShortAdaptation_IndvLeg(model{i},Ymuscles(:,:,i)',Uf,analysis,{labels(i).Data(2:end-1)},isF)
%     set(gcf,'color','w')
    
    
end
%% 
% %% Plotting all the muscles for individual participants
% 
% muscleOrder={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT',...
%     'VM', 'VL', 'RF', 'HIP','TFL', 'GLU'};
% ytl=([strcat('f',muscleOrder) strcat('s',muscleOrder)]);  %List of muscle
% ytl(end:-1:1) = ytl(:);
% yt=1:length(ytl);
% fs=14;
% c=colormap(jet(12));
% 
% %EMG Data
% 
% if type==1 %Testing
%     groupID='BATS'
% elseif type==2 %Training
%     groupID='BATR'
% elseif type==3 %Both groups
%     groupID='BAT'
% end
% 
% files = dir ([groupID '*params.mat']);
% 
% 
% n_subjects = size(files,1);
% 
% ii=0;
% for i =1:n_subjects
%     ii=1+ii;
%     sub{ii} = files(i).name;
%     subID{ii} = sub{ii}(1:end-10);
% end
% 
% figure()
% 
% for s=1:12
%     
%     subplot 211
%     hold on
%     
%     li{s}=scatter(1:28,reactive_trace(s,:),'filled','MarkerFaceColor', c(s,:));
%     
%     subplot 212
%     hold on
%     scatter(1:28,contextual_trace(s,:),'filled','MarkerFaceColor', c(s,:))
%     
% end
% 
% legend([li{:}],subID{:},'AutoUpdate','off')
% subplot 211
% hold on
% ylabel({'Reactive';'AU'})
% yline(0)
% set(gca,'XTick',yt,'XTickLabel',ytl,'FontSize',10)
% 
% subplot 212
% hold on
% set(gca,'XTick',yt,'XTickLabel',ytl,'FontSize',10)
% yline(0)
% ylabel({'Contextual';'AU'})
% 
% set(gcf,'color','w')
% 
% 
% %% Geting the average of the first 5 strides post-adaptation per muscle
% figure
% hold on
% for i=1:28
%     if i<15
%         Li{1}=scatter(nanmean(temp2(481:485,1,i)),nanmean(temp2(481:485,2,i)),100,"filled",'MarkerFaceColor', 'b');
%         text(nanmean(temp2(481:485,1,i)),nanmean(temp2(481:485,2,i)),{labels(i).Data(2:end-1)})
%     else
%         Li{2}=scatter(nanmean(temp2(481:485,1,i)),nanmean(temp2(481:485,2,i)),100,"filled",'MarkerFaceColor', 'r')  ;
%         text(nanmean(temp2(481:485,1,i)),nanmean(temp2(481:485,2,i)),{labels(i).Data(2:end-1)})
%     end
%     
% end
% legend([Li{:}],['Slow';'Fast'])
% ylabel({'Contextual';'A.U'})
% xlabel({'Reactive';'A.U'})
% 
% xlim([-0.3 1.4])
% ylim([-0.3 1.4])
% set(gcf,'color','w')
% 
% %% Geeting the average of Early Adapt and lateAdapt
% figure
% hold on
% for i=1:28
%     if i<15
%         Li{1}=scatter(nanmean(temp2(41:45,1,i)),nanmean(temp2(41:45,2,i)),100,"filled",'MarkerFaceColor', 'b');
%         text(nanmean(temp2(41:45,1,i)),nanmean(temp2(41:45,2,i)),{labels(i).Data(2:end-1)})
%     else
%         Li{2}=scatter(nanmean(temp2(41:45,1,i)),nanmean(temp2(41:45,2,i)),100,"filled",'MarkerFaceColor', 'r')  ;
%         text(nanmean(temp2(41:45,1,i)),nanmean(temp2(41:45,2,i)),{labels(i).Data(2:end-1)})
%     end
%     
% end
% legend([Li{:}],['Slow';'Fast'])
% ylabel({'Contextual';'A.U'})
% xlabel({'Reactive';'A.U'})
% % axis square
% xlim([-1 3])
% ylim([-1 1])
% 
% title('EarlyAdapt')
% set(gcf,'color','w')
% 
% 
% figure
% hold on
% for i=1:28
%     if i<15
%         Li{1}=scatter(nanmean(temp2(440:480,1,i)),nanmean(temp2(440:480,2,i)),100,"filled",'MarkerFaceColor', 'b');
%         text(nanmean(temp2(440:480,1,i)),nanmean(temp2(440:480,2,i)),{labels(i).Data(2:end-1)})
%     else
%         Li{2}=scatter(nanmean(temp2(440:480,1,i)),nanmean(temp2(440:480,2,i)),100,"filled",'MarkerFaceColor', 'r')  ;
%         text(nanmean(temp2(440:480,1,i)),nanmean(temp2(440:480,2,i)),{labels(i).Data(2:end-1)})
%     end
%     
% end
% legend([Li{:}],['Slow';'Fast'])
% ylabel({'Contextual';'A.U'})
% xlabel({'Reactive';'A.U'})
% axis square
% title('LateAdapt')
% set(gcf,'color','w')
% xlim([-.1 .85])
% ylim([-.1 .85])
% 
% %% Comparing individual analysis with the leg specific analysis
% 
% Cslow=[Casym(1:size(C,1)/2,1)  Casym(1:size(C,1)/2,2)]./vecnorm([Casym(1:size(C,1)/2,1)  Casym(1:size(C,1)/2,2)]); %SLOW
% Cfast=[Casym(1+size(C,1)/2:end,1)  Casym(1+size(C,1)/2:end,2)]./vecnorm([Casym(1+size(C,1)/2:end,1)  Casym(1+size(C,1)/2:end,2)]); %FAST
% 
% Ymodel=Ymodel';
% Yslow=Ymodel(1:680,1:size(Ymodel,2)/2); %SLOW
% Yfast=Ymodel(1:680,size(Ymodel,2)/2+1:end); %FAST
% 
% 
% 
% Cinv=pinv(Cslow);
% Xslow = Cinv*Yslow'; %x= y/C
% hatYslow=  Cslow *Xslow ; %yhat = C
% 
% % hatYslow= hatYslow';
% 
% Cinv=pinv(Cfast);
% Xfast = Cinv*Yfast'; %x= y/C
% hatYfast=  Cfast *Xfast ; %yhat = C
% 
% % hatYfast=hatYfast;
% Yfast=Yfast';
% Yslow=Yslow';
% 
% % Variance explained
% 
% ex2=[0.2314    0.2980    0.7529];
% ex1=[0.7255    0.0863    0.1608];
% mid=ones(1,3);
% N=100;
% gamma=1.5; %gamma > 1 expands the white (mid) part of the map, 'hiding' low values. Gamma<1 does the opposite
% gamma=1;
% map=[flipud(mid+ (ex1-mid).*([1:N]'/N).^gamma); mid; (mid+ (ex2-mid).*([1:N]'/N).^gamma)];
% 
% 
% 
% %% Variance all muscles
% Rsquared_slow2= my_Rsquared_coeff(Yslow,hatYslow);
% Rsquared_fast2 = my_Rsquared_coeff(Yfast,hatYfast);
% 
% Rsquared_slow_unit2 = my_Rsquared_coeff(Yslow,reconstruction_indv(1:168,:));
% Rsquared_fast_unit2 = my_Rsquared_coeff(Yfast,reconstruction_indv(169:end,:));
% 
% figure()
% plot(Rsquared_slow2)
% hold on
% plot(Rsquared_slow_unit2)
% legend('Per leg','Per muscle')
% 
% figure
% plot(Rsquared_fast2)
% hold on
% plot(Rsquared_fast_unit2)
% legend('Per leg','Per muscle')
% 
% % both legs
% Rsquared_both= my_Rsquared_coeff([Yslow(:,:);Yfast(:,:)],[hatYslow(:,:);hatYfast(:,:)]);
% Rsquared_both_unit = my_Rsquared_coeff([Yslow(:,:);Yfast(:,:)],reconstruction_indv(:,:));
% 
% figure
% plot(Rsquared_both)
% hold on
% plot(Rsquared_both_unit)
% legend('Per leg','Per muscle')
% title('Both leg')
% set(gcf,'color','w')
% % figure(3)
% % ytl={'GLU','TFL','HIP','RF','VL','VM','BF', 'SEMB','SEMT','MG','LG','SOL','PER','TA'};
% %%
% %% Calculating the R squared per muscle
% % mm= 0:12:168;
% % mm2=1:12:168;
% % for m=1:14
% %     Rsquared_slow(m,:) = my_Rsquared_coeff(Yslow(mm2(m):mm(m+1),:),hatYslow(mm2(m):mm(m+1),:));
% %     Rsquared_fast(m,:) = my_Rsquared_coeff(Yfast(mm2(m):mm(m+1),:),hatYfast(mm2(m):mm(m+1),:));
% %
% %     Rsquared_slow_unit(m,:) = my_Rsquared_coeff(Ymuscles(:,:,m)',Y2asym(:,:,m));
% %     Rsquared_fast_unit(m,:) = my_Rsquared_coeff(Ymuscles(:,:,m+14)',Y2asym(:,:,m+14));
% %
% % end
% % % figure
% % muscleOrder={'GLU','TFL','HIP','RF','VL','VM','BF', 'SEMB','SEMT','MG','LG','SOL','PER','TA'};
% % ytl= defineMuscleListV2(muscleOrder); %List of muscle
% % % ytl={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT', 'VM', 'VL', 'RF', 'HIP','TFL', 'GLU'};
% % % ytl(end:-1:1) = ytl(:);
% % binw=5;
% %
% % mtp=10;
% % for m=mtp
% %     %     figure
% %     %      imagesc(nanmean(Yslow(mm2(m):mm(m+1),481:485),2)',[-1 ,1])
% %     %      colormap(flipud(map))
% %     figure
% %     %     subplot(14,1,m)
% %     hold on
% %     aux1=conv(Rsquared_slow(m,:),ones(1,binw)/binw,'valid'); %Smoothing
% %     plot(aux1,'LineWidth',2,'DisplayName','Per leg','Color',"#0072BD") ;
% %
% %     aux1=conv(Rsquared_slow_unit(m,:),ones(1,binw)/binw,'valid'); %Smoothing
% %     plot(aux1,'LineWidth',2,'DisplayName','Per muscle','Color',"#A2142F") ;
% %
% %     %     plot(movmean(muscles_slow(m,:),5))
% %     %     plot(movmean(muscles_slow_unit(m,:),5))
% %     %     scatter(1:length(muscles(m,:)),movmean(muscles(m,:),5),'filled')
% %     %     legend({'Per leg';'Per muscle'},'AutoUpdate','off');
% %     legend('Location','NorthEastOutside','AutoUpdate','off')
% %     ylabel(ytl{m+14})
% %     %     yline(nanmean(muscles_slow(m,10:30)))
% %
% %     pp=patch([40 480 480 40],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
% %
% %     %     pp=patch([0 440  440 0],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
% %
% %     uistack(pp,'bottom')
% % end
% % set(gcf,'color','w')
% %
% % % figure
% % for m=mtp
% %     %     subplot(14,1,m)
% %     figure
% %     hold on
% %     aux1=conv(Rsquared_fast(m,:),ones(1,binw)/binw,'valid'); %Smoothing
% %     plot(aux1,'LineWidth',2,'DisplayName','Per leg','Color',"#0072BD") ;
% %
% %     aux1=conv(Rsquared_fast_unit(m,:),ones(1,binw)/binw,'valid'); %Smoothing
% %     plot(aux1,'LineWidth',2,'DisplayName','Per muscle','Color',"#A2142F") ;
% %     %         plot(movmean(muscles_fast(m,:),5))
% %     %         plot(movmean(muscles_fast_unit(m,:),5))
% %     %     scatter(1:length(muscles(m,:)),movmean(muscles(m,:),5),'filled')
% %     %     legend({'Per leg';'Per muscle'},'AutoUpdate','off');
% %
% %     legend('Location','NorthEastOutside','AutoUpdate','off')
% %     ylabel(ytl{m})
% %     %     yline(nanmean(muscles_fast(m,10:30)))
% %     pp=patch([40 480 480 40],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
% %     %     pp=patch([0 440  440 0],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
% %
% %     uistack(pp,'bottom')
% % end
% % set(gcf,'color','w')