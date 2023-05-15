% This is a script to get the confidance interval of the step-by-step
% weights

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
 

%% Load data and Plot checkerboard for all conditions.
% clear all; close all; clc;
% clear all; clc;

groupID ='BATS';
[normalizedGroupData, newLabelPrefix,n_subjects]=creatingGroupdataWnormalizedEMG(groupID);

%% Removing bad muscles
%This script make sure that we always remove the same muscle for the
%different analysis
normalizedGroupData= RemovingBadMuscleToSubj(normalizedGroupData);

%%  Getting the step-by-step data

% % Adaptation epochs
% strides=[-40 440 200];
%  
% if contains(groupID,'TR') %for treadmill Post 1
%     cond={'TM base','Adaptation','Post 1'}; %Conditions for this group
% else % for overground post 1
%     cond={'OG base','Adaptation','Post 1'}; %Conditions for this group
% end
% 
% exemptFirst=[1];
% exemptLast=[5]; %Strides needed
% names={};
% shortNames={};
% 
% ep=defineEpochs(cond,cond,strides,exemptFirst,exemptLast,'nanmean',{'Base','Adapt','Post1'}); %epochs
% 
% padWithNaNFlag=true; %If no enough steps fill with nan, let this on
% [dataEMG,labels,allDataEMG]=normalizedGroupData.getPrefixedEpochData(newLabelPrefix(end:-1:1),ep,padWithNaNFlag); %Getting the data
% 
% %Flipping EMG:
% for i=1:length(allDataEMG)
%     aux=reshape(allDataEMG{i},size(allDataEMG{i},1),size(labels,1),size(labels,2),size(allDataEMG{i},3));
%     allDataEMG{i}=reshape(flipEMGdata(aux,2,3),size(aux,1),numel(labels),size(aux,4));
% end
% 
% EMGdata=cell2mat(allDataEMG); %Getting EMG data per participants
% muscPhaseIdx=1:size(EMGdata,2); %


%% Getting the C values
% epochOfInterest={'TM base','TM mid 1','PosShort_{early}','PosShort_{late}','Ramp','Optimal','Adaptation','Adaptation_{early}','TiedPostPos','TMmid2','NegShort_{late}','Post1_{Early}','TMbase_{early}'};
epochOfInterest={'TM base','NegShort_{late}','Ramp','Optimal'};

ep=defineRegressorsDynamicsFeedback('nanmean');

if contains(groupID,'TR') %epoch to use for bias removal
    refEpTM = defineReferenceEpoch('TM base',ep);
else
    refEpTM = defineReferenceEpoch('OG base',ep);
end
flip=1;

if flip==1
    n=2;
    method='IndvLegs';
else
    n=1;
    method='Asym';
end

for s=1:n_subjects
    for l=1:length(epochOfInterest)
        ep2=defineReferenceEpoch(epochOfInterest{l},ep);
        adaptDataSubject = normalizedGroupData.adaptData{1, s};
        [~,~,~,Data{s,l}]=adaptDataSubject.getCheckerboardsData(newLabelPrefix,ep2,[],flip);
    end
end

%%

if strcmp(groupID,'BATS')
    fname='dynamicsData_BATS_subj_12_RemoveBadMuscles1_splits_0_WithPost2V2_WogBaseline.h5'
    load BATS_12_IndvLegsC17_ShortPertubations_RemovedBadMuscle_1RemovBias_0.mat
elseif  strcmp(groupID,'BATR')
    fname='dynamicsData_BATR_subj_12_RemoveBadMuscles1_splits_0_WithPost2V2.h5'
    load BATR_12_IndvLegsC16_ShortPertubations_RemovedBadMuscle_1RemovBias_0.mat
end

EMGdata2=h5read(fname,'/EMGdata');
 
binwith=10;
[~,~,~,~,~,~,EMGdata,labels]=groupDataToMatrixForm_Update(1:size(EMGdata2,3),fname,0);
muscPhaseIdx=1:size(EMGdata,2); %

context= find(strcmp(epochOfInterest,'Optimal')==1);
reactive2=find(strcmp(epochOfInterest,'NegShort_{late}')==1);



% mix=[C(1:168,context); C(169:end,reactive2)];
% Casym=[C(:,reactive2) C(:,context)]; % EMGreactive and EMGcontext
% Ymodel=Yasym';

%% Bootstrapping

epochOfInterest={'TM base','NegShort_{late}','Ramp','Optimal'};
context= find(strcmp(epochOfInterest,'Optimal')==1);
reactive2=find(strcmp(epochOfInterest,'NegShort_{late}')==1);

bootstrap=1; %Do you want to run the loop (1 yes 0 No)
X1=[];
X2=[];
replacement=1; %do you want to do it with replacement (1 yes 0 No)


if bootstrap
    if replacement
        n=2000; %number of iterations
    else
        n=1;
    end
    
    f = waitbar(0,'1','Name','Boostrapping Data',...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    
    setappdata(f,'canceling',0);
    
    unit=nan(12,2,28,n); %Creating nan matrices
    temp3=nan(12,2,28,n);
    Yhat=nan(12,680,28,n);
    temp4=nan(680,2,28,n);
    
    for l=1:n %loop for number of iterations
        
        temp=[];
        x=[];
        DataBoot={};
        % Check for clicked Cancel button
        if getappdata(f,'canceling')
            break
        end
        % Update waitbar and message
        ww=waitbar(l/n,f,['Iteration ' num2str(l)]);
        
        if replacement %doing the bootstrap with replacement
            subjIdx=datasample(1:n_subjects,n_subjects,'Replace',true);
        else
            subjIdx=datasample(1:n_subjects,n_subjects,'Replace',false);
        end
        
        DataBoot=Data(subjIdx,:); %Subject pick at each loop
        
        %This loop is to compute the constrant on our regressions
        for c=1:length(epochOfInterest)
            for s=1:n_subjects
                temp(:,:,s)=DataBoot{s,c};
            end
            x{1,c}=nanmedian(temp,3);
            tt=x{1,c}(:,end:-1:1);
            x{1,c}=tt;
            x{1,c}=reshape(x{1,c},14*2*12,1); %reshping the data for the C values 
        end
        
        x=cell2mat(x)';
        x=x';
        
        %C values that we are using for the regressions
        C2=[x(:,reactive2) x(:,context)];
        
        %Picking the data muscles that we want and participants
        Y2=EMGdata(:,muscPhaseIdx,subjIdx);
        
        %removing the bias for group
        Y2=nanmedian(Y2,3); %getting the median of the group 
        bias=nanmean(Y2(5:30,:,:)) ; %estimating the gorup baseline 
        C2=C2-bias'; %removing the bias from the constants 
        Y2=Y2-bias; %removing the bias from the data
        
        %reorganize the data to be separatend by muscle
        Cmuscles=reshape(C2',2,12,28);
        Ymuscles=reshape(Y2(1:680,:),680,12,28);
        
        %Linear regression individual muscles
        reconstruction_indv=[];
        data=[];
        C_indv=[];
        X_indv=[];
        
        for i=1:size(Ymuscles,3) %loop for individual muscle fit
            
            
            unit(:,:,i,l)=Cmuscles(:,:,i)'./vecnorm(Cmuscles(:,:,i)');
            %             unit(:,1,i,l)=-unit(:,1,i,l);
            temp3(:,:,i,l)=pinv(unit(:,:,i,l)'); %geeting the inverse of the constant
            Xhat(:,:,i,l) =temp3(:,:,i,l)'*Ymuscles(:,:,i)'; %x= y/C
            Yhat(:,:,i,l)=  unit(:,:,i,l)* Xhat(:,:,i,l) ; %Estimated Y with the constants 
            dynamics(:,:,i,l)=Xhat(:,:,i,l)'; %step-by-step dynamics
            
            
        end
        
    end
    close(ww)
    
end

delete(f)
%%
save([groupID,'_',num2str(n_subjects),'_iteration_', num2str(n),'_Individual_muscles'],'dynamics','Yhat','Ymuscles','groupID','-v7.3')
%%
load('musclesLabels.mat')
load('BATS_12_iteration_2000_Individual_muscles.mat')
OG=dynamics;
load('BATR_12_iteration_2000_Individual_muscles.mat')
TM=dynamics;


load('musclesLabels.mat')
load('BATR_indv_muscles.mat')
TM_2=X2asym;
load('BATS_indv_muscles.mat')
OG_2=X2asym;

%%
clrMap = colorcube(28*3);
muscles=[1:14 1:14];
g=[14:-1:1 14:-1:1];
ff=[1:14 1:14;15:28 15:28];
range=481:485;
for dyn=1:2%
    figure()
    hold on
    temp=[];
    x=[];
    
    for m=1:28
%         figure(ff(dyn,m))
%         hold on
        x=squeeze(TM(:,dyn,m,:));
        y=squeeze(OG(:,dyn,m,:));
        
        x_mean=nanmean(x(range,:),'all');
        y_mean=nanmean(y(range,:),'all');
        
        centers=[x_mean  y_mean];
        
        P_x = prctile(nanmean(x(range,:),1)',[2.5 97.5],"all");
        P_y = prctile(nanmean(y(range,:),1)',[2.5 97.5],"all");
        
        llc=[P_x(1), P_y(1)];
        
        CIrng(1)=P_x(2)-P_x(1);
        CIrng(2)=P_y(2)-P_y(1);
        
        %          CIrng=[CIrng_TM, CIrng_OG];
        %
        %         pd = fitdist(nanmean(temp(481:490,:),1)','Normal');
        %         ci = paramci(pd);
        %         pd_tm = fitdist(nanmean(temp2(481:490,:),1)','Normal');
        %         ci_tm = paramci(pd_tm);
%         OG_std=nanstd(y,0,'all');
%         TM_std=nanstd(x,0,'all');
        
        %             scatter(OG_mean,TM_mean,'fillled')
        %                 errorbar(TM_mean,OG_mean,OG_mean-P_OG(1),OG_mean-P_OG(2),TM_mean-P_TM(1),TM_mean-P_TM(2),"o")
        
        %         a= ci_tm(2) ; %CIrng_TM; % horizontal radius
        %         b=ci(2) ; %CIrng_OG; % vertical radius
        x0=x_mean; % x0,y0 ellipse centre coordinates
        y0=y_mean;
        %         t=-pi:0.01:pi;
        %         x=x0+a*cos(t);
        %         y=y0+b*sin(t);
        %         plot(x,y)
        
        text(x0+.02,y0,{labels(m).Data(1:end-1)})
        if P_y(1)<0 &&  P_y(2)>0 %P_TM(1)<0 &&  P_TM(2)>0 ||
            rectangle('Position',[llc,CIrng],'Curvature',[1,1],'EdgeColor',clrMap(m+3,:),'LineStyle','--');
        else
            rectangle('Position',[llc,CIrng],'Curvature',[1,1],'EdgeColor',clrMap(m+3,:));
        end

        
        plot(centers(1), centers(2), 'o', 'MarkerFaceColor', clrMap(m+3,:), 'MarkerSize',10, 'LineWidth', 1,'MarkerEdgeColor',clrMap(m+3,:))%clrMap(m+3,:))
        xlabel('Treadmill')
        ylabel('Overground')
%         
%         if m<15
%             Li{1}=scatter(nanmean(TM_2(dyn,range,m)),nanmean(OG_2(dyn,range,m)),100,"filled",'MarkerFaceColor', 'b');
%             text(nanmean(TM_2(dyn,range,m))+.02,nanmean(OG_2(dyn,range,m)),{labels(m).Data(1:end-1)})
%         else
%             Li{2}=scatter(nanmean(TM_2(dyn,range,m)),nanmean(OG_2(dyn,range,m)),100,"filled",'MarkerFaceColor', 'r')  ;
%                text(nanmean(TM_2(dyn,range,m))+.02,nanmean(OG_2(dyn,range,m)),{labels(m).Data(1:end-1)})
%         end

% if dyn==1
%     title('Reactive')
%     xx=-.5:0.1:2.1;
%     plot(xx,xx,'r')
%     xlim([-.5 2.1])
%     ylim([-.5 2.1])
%     yline(0)
%     xline(0)
% else
%     title('Contextual')
%     yline(0)
%     xline(0)
%     xlim([-.5 1])
%     ylim([-.5 1])
% end
% xlabel('Treadmill')
% ylabel('Overground')
%     set(gcf,'color','w')
    end
    xlabel('Treadmill')
    ylabel('Overground')
%     
    if dyn==1
        title('Reactive')
        xx=-.5:0.1:2.1;
        plot(xx,xx,'r')
        xlim([-.5 2.1])
        ylim([-.5 2.1])
        yline(0)
        xline(0)
    else
        title('Contextual')
        yline(0)
        xline(0)
        xlim([-.5 1.5])
        ylim([-.5 1.5])
    end
    
    set(gcf,'color','w')
    %     figure(dyn+10)
    %     modifiedBoxPlot([1:2],[nanmean(x(481:490,:),1)' nanmean(y(481:490,:),1)'])
    %     set(gca,'XTick',[1 2],'XTickLabel',{'Reactive^{OG}_{earlyPost}','Reactive^{TM}_{earlyPost}' },'FontSize',10)
    
end



%%
grayColor = [.7 .7 .7];

figure()
subplot(2,1,1)
hold on
% X_1=[X1(:,1:40) nan(size(X1(:,1:40),1),1) X1(:,41:480) nan(size(X1(:,1:40),1),1) X1(:,481:end)];
% y = nanmean(X1,1); % your mean vector;
% x = 1:numel(y);
% std_dev = nanstd(X1,1);
% curve1 = y + std_dev;
% curve2 = y - std_dev;
% x2 = [x, fliplr(x)];
% inBetween = [curve1, fliplr(curve2)];
% fill(x2, inBetween,'r','FaceAlpha',0.3,'EdgeColor','none')
% hold on;
% ylabel('W_{reactive}')
% xlabel('Strides')
% plot(x, y, 'r', 'LineWidth', 2);

% d{1} = X1(:,1:40);
% d{2} = X1(:,41:480);
% d{3} =  X1(:,481:end);
% index{1}= 1:40;
% index{2}= 42:481;
% index{3}= 483:682;
%
% for i=1:3
%     y = nanmean(d{i},1); % your mean vector;
%     %     x = 1:numel(y);
%     x=index{i};
%     std_dev = nanstd(d{i},1);
%     curve1 = y + std_dev;
%     curve2 = y - std_dev;
%     x2 = [x, fliplr(x)];
%     inBetween = [curve1, fliplr(curve2)];
%     fill(x2, inBetween,'r','FaceAlpha',0.3,'EdgeColor','none')
%     hold on;
%     ylabel('W_{reactive}')
%     xlabel('Strides')
%     plot(x, y, 'r', 'LineWidth', 2);
%
% end
% yline(0)
% % legend('Training')
%
%
% subplot(2,1,2)
% hold on
% % X_2=[X2(:,1:40) nan(size(X2(:,1:40),1),1) X2(:,41:480) nan(size(X2(:,1:40),1),1) X2(:,481:end)];
% % y = nanmean(X2,1); % your mean vector;
% % x = 1:numel(y);
% % std_dev = nanstd(X2,1);
% % curve1 = y + std_dev;
% % curve2 = y - std_dev;
% % x2 = [x, fliplr(x)];
% % inBetween = [curve1, fliplr(curve2)];
% % fill(x2, inBetween,'r','FaceAlpha',.3,'EdgeColor','none')
% % hold on;
% % ylabel('W_{context}')
% % xlabel('Strides')
% % plot(x, y, 'r', 'LineWidth', 2);
% d{1} = X2(:,1:40);
% d{2} = X2(:,41:480);
% d{3} =  X2(:,481:end);
%
% for i=1:3
%     y = nanmean(d{i},1); % your mean vector;
%     %     x = 1:numel(y);
%     x=index{i};
%     std_dev = nanstd(d{i},1);
%     curve1 = y + std_dev;
%     curve2 = y - std_dev;
%     x2 = [x, fliplr(x)];
%     inBetween = [curve1, fliplr(curve2)];
%     fill(x2, inBetween,'r','FaceAlpha',0.3,'EdgeColor','none')
%     hold on;
%     ylabel('W_{context}')
%     xlabel('Strides')
%     plot(x, y, 'r', 'LineWidth', 2);
%
% end
% yline(0)
% set(gcf,'color','w')
% %%
% save([groupID,'_',num2str(n_subjects),'_iteration_', num2str(n)],'dynamics','Yhat','Ymuscles','groupID')
%
%
% %% EMG norm
% load('BATR_12_iteration_10000wEMGnorm_16-Feb-2023.mat')
% subplot(2,1,1)
% hold on
% TR_post=nanmean(EMGnorm(:,481:485),2);
% TR_eA=nanmean(EMGnorm(:,41:51),2);
% TR_lA=nanmean(EMGnorm(:,480-40:480),2);
% histogram(TR_post,'FaceColor','r')
% xlabel('EMGnorm')
% disp('BATR')
% % [h,p,ci,stats] = ttest(reactive)
% % subplot(2,1,2)
% hold on
%
% load('BATS_12_iteration_10000wEMGnorm_16-Feb-2023.mat')
% TS_post=nanmean(EMGnorm(:,481:485),2);
% TS_eA=nanmean(EMGnorm(:,41:51),2);
% TS_lA=nanmean(EMGnorm(:,480-40:480),2);
% histogram(TS_post,'FaceColor','b')
% disp('BATR - Context')
% % [h,p,ci,stats]= ttest(context)
% xlabel('Context')
% legend('BATR','BATS')
% set(gcf,'color','w')
%
% disp('Context Distribution test')
% [h,p,ks2stat] = kstest2( TS_post, TR_post)
%
% % disp('t-test')
% % [h,p,ks2stat] = ttest2( TS_post, TR_post)
%
% title(['Two-sample Kolmogorov-Smirnov, p=',num2str(p),' Rejects the null hypothesis =', num2str(h)])
%
% stats_norm=modifiedBoxPlot([1,2],[ TS_post TR_post]);
% % set(gca,'title','Marce')
% set(gca,'XTick',[1 2],'XTickLabel',{'EMGnorm^{TS}_{post}','EMGnorm^{TR}_{post}' },'FontSize',10)
% title('Post-adaptation')
%
% stats_norm=modifiedBoxPlot([1,2,3,4],[ TS_eA TR_eA TS_lA TR_lA]);
% % set(gca,'title','Marce')
% set(gca,'XTick',[1 2,3,4],'XTickLabel',{'EMGnorm^{TS}_{early}','EMGnorm^{TR}_{early}' ,'EMGnorm^{TS}_{late}','EMGnorm^{TR}_{late}' },'FontSize',10)
% title('Adaptation')
%
% %% Adaptation
%
% eA=42:51;
% lA=435:480-5;
% eP=482:486;
%
%
% % load('BAT_24_iteration_10000wEMGnorm_16-Feb-2023.mat')
% % figure
% % subplot(2,1,1)
% % hold on
% % early_reactive1=nanmean(X1(:,eA),2);
% % late_reactive1=nanmean(X1(:,lA),2);
% %
% %
% % early_context1=nanmean(X2(:,eA),2);
% % late_context1=nanmean(X2(:,lA),2);
% % histogram(early_context1)
% %
% %
% % stats_Adapt=modifiedBoxPlot([1:4],[early_reactive1  early_context1   late_reactive1 late_context1 ]);
% % set(gca,'XTick',[1 2,3,4],'XTickLabel',{'W_{early-reactive}','W_{early-contex}','W_{late-reactive}','W_{late-contex}'},'FontSize',10)
% % title('Adaptation (All Data)')
% %
% % %% Post-adaptation data
% % eA=42:51;
% % lA=435:480-5;
% % eP=482:486;
% % lP=680-45:680-5;
% %
% % % load('BATS_12_iteration_10000wEMGnorm_16-Feb-2023.mat')
% % load('BATS_12_iteration_10000wEMGnorm_TiedPostRamp13-Mar-2023')
% % figure
% % subplot(2,1,1)
% % hold on
% %
% % early_reactive1=nanmean(X1(:,eA),2);
% % late_reactive1=nanmean(X1(:,lA),2);
% % reactive1=nanmean(X1(:,eP),2);
% % late_Preactive1=nanmean(X1(:,lP),2);
% % histogram(reactive1)
% % xlabel('Reactive')
% % disp('BATS - Reactive')
% %
% % [h,p,ci,stats] = ttest(reactive1)
% %
% % subplot(2,1,2)
% % hold on
% % context1=nanmean(X2(:,eP),2);
% % early_context1=nanmean(X2(:,eA),2);
% % late_context1=nanmean(X2(:,lA),2);
% % late_Pcontext1=nanmean(X2(:,lP),2);
% % histogram(context1)
% % disp('BATS - Context')
% % [h,p,ci,stats]= ttest(context1)
% % xlabel('Context')
% %
% %
% % % early_reactive2=nanmean(X3(:,eA),2);
% % % late_reactive2=nanmean(X3(:,lA),2);
% % % reactive2=nanmean(X3(:,eP),2);
% % % late_Preactive2=nanmean(X3(:,lP),2);
% %
% % load('BATR_12_iteration_10000wEMGnorm_TiedPostRamp13-Mar-2023')
% % % load('BATR_12_iteration_10000wEMGnorm_16-Feb-2023.mat')
% % % load('BATR_11_Individual_C.mat')
% % figure
% % subplot(2,1,1)
% % hold on
% % reactive=nanmean(X1(:,eP),2);
% % early_reactive=nanmean(X1(:,eA),2);
% % late_reactive=nanmean(X1(:,lA),2);
% % late_Preactive=nanmean(X1(:,lP),2);
% % %
% % % reactive3=nanmean(X3(:,eP),2);
% % % early_reactive3=nanmean(X3(:,eA),2);
% % % late_reactive3=nanmean(X3(:,lA),2);
% % % late_Preactive3=nanmean(X3(:,lP),2);
% %
% % histogram(reactive)
% % xlabel('Reactive')
% % disp('BATR - Reactive')
% % [h,p,ci,stats] = ttest(reactive)
% % subplot(2,1,2)
% % hold on
% % context=nanmean(X2(:,eP),2);
% % early_context=nanmean(X2(:,eA),2);
% % late_context=nanmean(X2(:,lA),2);
% % late_Pcontext=nanmean(X2(:,lP),2);
% % histogram(context)
% % disp('BATR - Context')
% % [h,p,ci,stats]= ttest(context)
% % xlabel('Context')
% % legend('BATS','BATR')
% % set(gcf,'color','w')
% %
% %
% % figure
% % % set(gca,'TickLabelInterpreter','latex');
% %
% % % stats_eP=modifiedBoxPlot([1:6],[reactive1 reactive reactive2 reactive3 context1 context]);
% % % set(gca,'XTick',[1 2,3,4,5,6],'XTickLabel',{'W^{OG}_{intro}','W^{TM}_{intro}','W^{OG}_{removal}','W^{TM}_{removal}','W^{OG}_{contex}','W^{TM}_{contex}'},'FontSize',10)
% % stats_eP=modifiedBoxPlot([1:4],[-reactive1 -reactive context1 context]);
% % set(gca,'XTick',[1 2,3,4],'XTickLabel',{'W^{OG}_{intro}','W^{TM}_{intro}','W^{OG}_{contex}','W^{TM}_{contex}'},'FontSize',10)
% %
% % title('Early Post-adaptation')
% % stats_eP.actualCI(1,:)=stats_eP.mean-abs(stats_eP.CI(1,:));
% % stats_eP.actualCI(2,:)=stats_eP.mean+abs(stats_eP.CI(2,:));
% %
% % display('Reactive Distribution test')
% % [h,p,ks2stat] = kstest2(reactive1,reactive)
% %
% %
% % display('Context Distribution test')
% % [h,p,ks2stat] = kstest2(context1, context)
% %
% %
% % % stats_lP=modifiedBoxPlot([1:6],[late_Preactive1 late_Preactive late_Preactive2 late_Preactive3 late_Pcontext1 late_Pcontext]);
% % % set(gca,'XTick',[1 2,3,4,5,6],'XTickLabel',{'W^{OG}_{intro}','W^{TM}_{intro}','W^{OG}_{removal}','W^{TM}_{removal}','W^{TM}_{contex}'},'FontSize',10)
% % stats_lP=modifiedBoxPlot([1:4],[-late_Preactive1 -late_Preactive late_Pcontext1 late_Pcontext]);
% % set(gca,'XTick',[1 2,3,4,],'XTickLabel',{'W^{OG}_{intro}','W^{TM}_{intro}','W^{TM}_{contex}'},'FontSize',10)
% %
% % title('Late Post-adaptation')
% %
% % stats_lP.actualCI(1,:)=stats_lP.mean-abs(stats_lP.CI(1,:));
% % stats_lP.actualCI(2,:)=stats_lP.mean+abs(stats_lP.CI(2,:));
% %
% % %% Expert behavior
% %
% % lE=1030-45:1030-5;
% %
% % % load('BATS_12_iteration_10000wEMGnorm_16-Feb-2023.mat')
% % figure% This is a script to get the confidance interval of the step-by-step
% % weights
% 
% %% Load data and Plot checkerboard for all conditions.
% % clear all; close all; clc;
% clear all; clc;
% 
% groupID ='BATS';
% [normalizedGroupData, newLabelPrefix,n_subjects]=creatingGroupdataWnormalizedEMG(groupID);
% 
% %% Removing bad muscles
% %This script make sure that we always remove the same muscle for the
% %different analysis
% normalizedGroupData= RemovingBadMuscleToSubj(normalizedGroupData);
% 
% %%  Getting the step-by-step data
% 
% %Adaptation epochs
% strides=[-40 440 200];
% 
% if contains(groupID,'TR') %for treadmill Post 1
%     cond={'TM base','Adaptation','Post 1'}; %Conditions for this group
% else % for overground post 1
%     cond={'OG base','Adaptation','Post 1'}; %Conditions for this group
% end
% 
% exemptFirst=[1];
% exemptLast=[5]; %Strides needed
% names={};
% shortNames={};
% 
% ep=defineEpochs(cond,cond,strides,exemptFirst,exemptLast,'nanmean',{'Base','Adapt','Post1'}); %epochs
% 
% padWithNaNFlag=true; %If no enough steps fill with nan, let this on
% [dataEMG,labels,allDataEMG]=normalizedGroupData.getPrefixedEpochData(newLabelPrefix(end:-1:1),ep,padWithNaNFlag); %Getting the data
% 
% %Flipping EMG:
% for i=1:length(allDataEMG)
%     aux=reshape(allDataEMG{i},size(allDataEMG{i},1),size(labels,1),size(labels,2),size(allDataEMG{i},3));
%     allDataEMG{i}=reshape(flipEMGdata(aux,2,3),size(aux,1),numel(labels),size(aux,4));
% end
% 
% EMGdata=cell2mat(allDataEMG); %Getting EMG data per participants
% muscPhaseIdx=1:size(EMGdata,2); %
% 
% 
% %% Getting the C values
% % epochOfInterest={'TM base','TM mid 1','PosShort_{early}','PosShort_{late}','Ramp','Optimal','Adaptation','Adaptation_{early}','TiedPostPos','TMmid2','NegShort_{late}','Post1_{Early}','TMbase_{early}'};
% epochOfInterest={'TM base','NegShort_{late}','Ramp','Optimal'};
% 
% ep=defineRegressorsDynamicsFeedback('nanmean');
% 
% if contains(groupID,'TR') %epoch to use for bias removal
%     refEpTM = defineReferenceEpoch('TM base',ep);
% else
%     refEpTM = defineReferenceEpoch('OG base',ep);
% end
% flip=1;
% 
% if flip==1
%     n=2;
%     method='IndvLegs';
% else
%     n=1;
%     method='Asym';
% end
% 
% for s=1:n_subjects
%     for l=1:length(epochOfInterest)
%         ep2=defineReferenceEpoch(epochOfInterest{l},ep);
%         adaptDataSubject = normalizedGroupData.adaptData{1, s};
%         [~,~,~,Data{s,l}]=adaptDataSubject.getCheckerboardsData(newLabelPrefix,ep2,refEpTM,flip);
%     end
% end
% 
% %% Bootstrapping
% 
% bootstrap=1; %Do you want to run the loop (1 yes 0 No)
% X1=[];
% X2=[];
% replacement=1; %do you want to do it with replacement (1 yes 0 No)
% 
% 
% if bootstrap
%     if replacement
%         n=2000; %number of iterations
%     else
%         n=3;
%     end
%     
%     f = waitbar(0,'1','Name','Boostrapping Data',...
%         'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
%     
%     setappdata(f,'canceling',0);
%     
%     unit=nan(12,2,28,n); %Creating nan matrices
%     temp3=nan(12,2,28,n);
%     Yhat=nan(12,680,28,n);
%     temp4=nan(680,2,28,n);
%     
%     for l=1:n %loop for number of iterations
%         
%         temp=[];
%         x=[];
%         DataBoot={};
%         % Check for clicked Cancel button
%         if getappdata(f,'canceling')
%             break
%         end
%         % Update waitbar and message
%         ww=waitbar(l/n,f,['Iteration ' num2str(l)]);
%         
%         if replacement %doing the bootstrap with replacement
%             subjIdx=datasample(1:n_subjects,n_subjects,'Replace',true);
%         else
%             subjIdx=datasample(1:n_subjects,n_subjects,'Replace',false);
%         end
%         
%         DataBoot=Data(subjIdx,:); %Subject pick at each loop
%         
%         for c=1:length(epochOfInterest)
%             for s=1:n_subjects
%                 temp(:,:,s)=DataBoot{s,c};
%             end
%             x{1,c}=nanmedian(temp,3);
%             tt=x{1,c}(:,end:-1:1);
%             x{1,c}=tt;
%             x{1,c}=reshape(x{1,c},14*2*12,1);
%         end
%         
%         x=cell2mat(x)';
%         x=x';
%         
%         %C values that we are using for the regressions
%         C2=[x(:,2) x(:,4)];
%         
%         %Picking the data muscles that we want and participants
%         Y2=EMGdata(:,muscPhaseIdx,subjIdx);
%         
%         %removing the bias for group
%         bias=nanmean(Y2(5:30,:,:)) ;
%         bias=nanmean(bias,3);
%         Y2=nanmedian(Y2,3)- bias;
%         
%         %reorganize the data to be separatend by muscle
%         Cmuscles=reshape(C2',2,12,28);
%         Ymuscles=reshape(Y2(1:680,:),680,12,28);
%         
%         %Linear regression individual muscles
%         reconstruction_indv=[];
%         data=[];
%         C_indv=[];
%         X_indv=[];
%         
%         for i=1:size(x,1)/12 %loop for individual muscle fit
%             
%             
%             unit(:,:,i,l)=Cmuscles(:,:,i)'./vecnorm(Cmuscles(:,:,i)');
%             %             unit(:,1,i,l)=-unit(:,1,i,l);
%             %         temp(:,:,i)=pinv(Cmuscles(:,:,i)');
%             temp3(:,:,i,l)=pinv(unit(:,:,i,l)');
%             Xhat(:,:,i,l) =temp3(:,:,i,l)'*Ymuscles(:,:,i)'; %x= y/C
%             Yhat(:,:,i,l)=  unit(:,:,i,l)* Xhat(:,:,i,l) ;
%             dynamics(:,:,i,l)=Xhat(:,:,i,l)';
%             
%             %            if mod(i,2)==0
%             %            figure(i)
%             %
%             %            subplot(2,1,1)
%             %             hold on
%             %            plot(Xhat(1,481:end,i,l))
%             %
%             %
%             %            subplot(2,1,2)
%             %            hold on
%             %            plot(Xhat(2,481:end,i,l))
%             %            end
%             
%             %             model{i}.C=unit(:,:,i,l);
%             
%             %             reconstruction_indv =[ reconstruction_indv ; Y2asym(:,:,i)];
%             %             data =[ data ; Ymuscles(:,:,i)'];
%             %             C_indv=[C_indv;unit(:,:,i)];
%             %             X_indv=[X_indv,X2asym(:,:,i) ];
%             %             %         model{i}.C=Cmuscles(:,:,i)';
%             %             temp5=corrcoef(model{i}.C);
%             %             correlation(i,1)=temp5(2);
%             %             temp5=vif([model{i}.C]);
%             %             impact(i,1)=temp5(1);
%             
%             
%         end
%         %         Yasym2=Y2-fftshift(Y2,2);
%         %         Yasym2=Yasym2(:,1:size(Yasym2,2)/2,:);
%         
%         
%         
%         %
%         %         Cinv2=pinv(C2);
%         %         Xdynamics= Cinv2*Y2'; %x= y/C
%         %
%         %         X1=[X1;Xdynamics(1,:)];
%         %         X2=[X2;Xdynamics(2,:)];
%         
%         
%         
%         
%     end
%     close(ww)
%     
% end
% 
% delete(f)
% %%
% save([groupID,'_',num2str(n_subjects),'_iteration_', num2str(n),'_Individual_muscles'],'dynamics','Yhat','Ymuscles','groupID','-v7.3')
% %%
% load('musclesLabels.mat')
% load('BATS_12_iteration_2000_Individual_muscles.mat')
% OG=dynamics;
% load('BATR_12_iteration_2000_Individual_muscles.mat')
% TM=dynamics;
% c                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
% 
% load('musclesLabels.mat')
% load('BATR_indv_muscles.mat')
% TM_2=X2asym;
% load('BATS_indv_muscles.mat')
% OG_2=X2asym;
% 
% %%
% clrMap = colorcube(28*3);
% 
% for dyn=1:2
% %     figure()
% %     hold on
%     temp=[];
%     x=[];
%     
%     for m=1:28
%         
%         figure()
%         hold on
%         %            figure()
%         %     hold on
%         x=squeeze(TM(:,dyn,m,:));
%         y=squeeze(OG(:,dyn,m,:));
%         
%         x_mean=nanmedian(x(481:490,:),'all');
%         y_mean=nanmedian(y(481:490,:),'all');
%         
%         centers=[x_mean  y_mean];
%         
%         P_x = prctile(nanmean(x(481:490,:),1)',[2.5 97.5],"all");
%         P_y = prctile(nanmean(y(481:490,:),1)',[2.5 97.5],"all");
%         
%         llc=[P_x(1), P_y(1)];
%         
%         CIrng(1)=P_x(2)-P_x(1);
%         CIrng(2)=P_y(2)-P_y(1);
%         
%         %          CIrng=[CIrng_TM, CIrng_OG];
%         %
%         %         pd = fitdist(nanmean(temp(481:490,:),1)','Normal');
%         %         ci = paramci(pd);
%         %         pd_tm = fitdist(nanmean(temp2(481:490,:),1)','Normal');
%         %         ci_tm = paramci(pd_tm);
%         OG_std=nanstd(y,0,'all');
%         TM_std=nanstd(x,0,'all');
%         
%         %             scatter(OG_mean,TM_mean,'fillled')
%         %                 errorbar(TM_mean,OG_mean,OG_mean-P_OG(1),OG_mean-P_OG(2),TM_mean-P_TM(1),TM_mean-P_TM(2),"o")
%         
%         %         a= ci_tm(2) ; %CIrng_TM; % horizontal radius
%         %         b=ci(2) ; %CIrng_OG; % vertical radius
%         x0=x_mean; % x0,y0 ellipse centre coordinates
%         y0=y_mean;
%         %         t=-pi:0.01:pi;
%         %         x=x0+a*cos(t);
%         %         y=y0+b*sin(t);
%         %         plot(x,y)
%         text(x0+.02,y0,{labels(m).Data(1:end-1)})
%         if P_y(1)<0 &&  P_y(2)>0 %P_TM(1)<0 &&  P_TM(2)>0 ||
%             rectangle('Position',[llc,CIrng],'Curvature',[1,1],'EdgeColor',clrMap(m+3,:),'LineStyle','--');
%         else
%             rectangle('Position',[llc,CIrng],'Curvature',[1,1],'EdgeColor',clrMap(m+3,:));
%         end
%         %         scatter(nanmean(temp2(481:490,:),1)',nanmean(temp(481:490,:),1)','k')
%         
%         plot(centers(1), centers(2), 'o', 'MarkerFaceColor', clrMap(m+3,:), 'MarkerSize',10, 'LineWidth', 1,'MarkerEdgeColor',clrMap(m+3,:))
%         xlabel('Treadmill')
%         ylabel('Overground')
%         
%         scatter(nanmean(TM_2(1,481:485,i)),nanmean(OG_2(1,481:485,i)),100,"filled",'MarkerFaceColor', 'b');
%         
%     end
%     xlabel('Treadmill')
%     ylabel('Overground')
%     
%     if dyn==1
%         title('Reactive')
%         xx=-.5:0.1:2.1;
%         plot(xx,xx,'r')
%         xlim([-.5 2.1])
%         ylim([-.5 2.1])
%         yline(0)
%         xline(0)
%     else
%         title('Contextual')
%         yline(0)
%         xline(0)
%         xlim([-.5 1.5])
%         ylim([-.5 1.5])
%     end
%     
%     set(gcf,'color','w')
%     %     figure(dyn+10)
%     %     modifiedBoxPlot([1:2],[nanmean(x(481:490,:),1)' nanmean(y(481:490,:),1)'])
%     %     set(gca,'XTick',[1 2],'XTickLabel',{'Reactive^{OG}_{earlyPost}','Reactive^{TM}_{earlyPost}' },'FontSize',10)
%     
% end
% 
% 
% 
% %%
% grayColor = [.7 .7 .7];
% 
% figure()
% subplot(2,1,1)
% hold on
% X_1=[X1(:,1:40) nan(size(X1(:,1:40),1),1) X1(:,41:480) nan(size(X1(:,1:40),1),1) X1(:,481:end)];
% y = nanmean(X1,1); % your mean vector;
% x = 1:numel(y);
% std_dev = nanstd(X1,1);
% curve1 = y + std_dev;
% curve2 = y - std_dev;
% x2 = [x, fliplr(x)];
% inBetween = [curve1, fliplr(curve2)];
% fill(x2, inBetween,'r','FaceAlpha',0.3,'EdgeColor','none')
% hold on;
% ylabel('W_{reactive}')
% xlabel('Strides')
% plot(x, y, 'r', 'LineWidth', 2);

% d{1} = X1(:,1:40);
% d{2} = X1(:,41:480);
% d{3} =  X1(:,481:end);
% index{1}= 1:40;
% index{2}= 42:481;
% index{3}= 483:682;
%
% for i=1:3
%     y = nanmean(d{i},1); % your mean vector;
%     %     x = 1:numel(y);
%     x=index{i};
%     std_dev = nanstd(d{i},1);
%     curve1 = y + std_dev;
%     curve2 = y - std_dev;
%     x2 = [x, fliplr(x)];
%     inBetween = [curve1, fliplr(curve2)];
%     fill(x2, inBetween,'r','FaceAlpha',0.3,'EdgeColor','none')
%     hold on;
%     ylabel('W_{reactive}')
%     xlabel('Strides')
%     plot(x, y, 'r', 'LineWidth', 2);
%
% end
% yline(0)
% % legend('Training')
%
%
% subplot(2,1,2)
% hold on
% % X_2=[X2(:,1:40) nan(size(X2(:,1:40),1),1) X2(:,41:480) nan(size(X2(:,1:40),1),1) X2(:,481:end)];
% % y = nanmean(X2,1); % your mean vector;
% % x = 1:numel(y);
% % std_dev = nanstd(X2,1);
% % curve1 = y + std_dev;
% % curve2 = y - std_dev;
% % x2 = [x, fliplr(x)];
% % inBetween = [curve1, fliplr(curve2)];
% % fill(x2, inBetween,'r','FaceAlpha',.3,'EdgeColor','none')
% % hold on;
% % ylabel('W_{context}')
% % xlabel('Strides')
% % plot(x, y, 'r', 'LineWidth', 2);
% d{1} = X2(:,1:40);
% d{2} = X2(:,41:480);
% d{3} =  X2(:,481:end);
%
% for i=1:3
%     y = nanmean(d{i},1); % your mean vector;
%     %     x = 1:numel(y);
%     x=index{i};
%     std_dev = nanstd(d{i},1);
%     curve1 = y + std_dev;
%     curve2 = y - std_dev;
%     x2 = [x, fliplr(x)];
%     inBetween = [curve1, fliplr(curve2)];
%     fill(x2, inBetween,'r','FaceAlpha',0.3,'EdgeColor','none')
%     hold on;
%     ylabel('W_{context}')
%     xlabel('Strides')
%     plot(x, y, 'r', 'LineWidth', 2);
%
% end
% yline(0)
% set(gcf,'color','w')
% %%
% save([groupID,'_',num2str(n_subjects),'_iteration_', num2str(n)],'dynamics','Yhat','Ymuscles','groupID')
%
%
% %% EMG norm
% load('BATR_12_iteration_10000wEMGnorm_16-Feb-2023.mat')
% subplot(2,1,1)
% hold on
% TR_post=nanmean(EMGnorm(:,481:485),2);
% TR_eA=nanmean(EMGnorm(:,41:51),2);
% TR_lA=nanmean(EMGnorm(:,480-40:480),2);
% histogram(TR_post,'FaceColor','r')
% xlabel('EMGnorm')
% disp('BATR')
% % [h,p,ci,stats] = ttest(reactive)
% % subplot(2,1,2)
% hold on
%
% load('BATS_12_iteration_10000wEMGnorm_16-Feb-2023.mat')
% TS_post=nanmean(EMGnorm(:,481:485),2);
% TS_eA=nanmean(EMGnorm(:,41:51),2);
% TS_lA=nanmean(EMGnorm(:,480-40:480),2);
% histogram(TS_post,'FaceColor','b')
% disp('BATR - Context')
% % [h,p,ci,stats]= ttest(context)
% xlabel('Context')
% legend('BATR','BATS')
% set(gcf,'color','w')
%
% disp('Context Distribution test')
% [h,p,ks2stat] = kstest2( TS_post, TR_post)
%
% % disp('t-test')
% % [h,p,ks2stat] = ttest2( TS_post, TR_post)
%
% title(['Two-sample Kolmogorov-Smirnov, p=',num2str(p),' Rejects the null hypothesis =', num2str(h)])
%
% stats_norm=modifiedBoxPlot([1,2],[ TS_post TR_post]);
% % set(gca,'title','Marce')
% set(gca,'XTick',[1 2],'XTickLabel',{'EMGnorm^{TS}_{post}','EMGnorm^{TR}_{post}' },'FontSize',10)
% title('Post-adaptation')
%
% stats_norm=modifiedBoxPlot([1,2,3,4],[ TS_eA TR_eA TS_lA TR_lA]);
% % set(gca,'title','Marce')
% set(gca,'XTick',[1 2,3,4],'XTickLabel',{'EMGnorm^{TS}_{early}','EMGnorm^{TR}_{early}' ,'EMGnorm^{TS}_{late}','EMGnorm^{TR}_{late}' },'FontSize',10)
% title('Adaptation')
%
% %% Adaptation
%
% eA=42:51;
% lA=435:480-5;
% eP=482:486;
%
%
% load('BAT_24_iteration_10000wEMGnorm_16-Feb-2023.mat')
% figure
% subplot(2,1,1)
% hold on
% early_reactive1=nanmean(X1(:,eA),2);
% late_reactive1=nanmean(X1(:,lA),2);
%
%
% early_context1=nanmean(X2(:,eA),2);
% late_context1=nanmean(X2(:,lA),2);
% histogram(early_context1)
%
%
% stats_Adapt=modifiedBoxPlot([1:4],[early_reactive1  early_context1   late_reactive1 late_context1 ]);
% set(gca,'XTick',[1 2,3,4],'XTickLabel',{'W_{early-reactive}','W_{early-contex}','W_{late-reactive}','W_{late-contex}'},'FontSize',10)
% title('Adaptation (All Data)')
%
% %% Post-adaptation data
% eA=42:51;
% lA=435:480-5;
% eP=482:486;
% lP=680-45:680-5;
%
% % load('BATS_12_iteration_10000wEMGnorm_16-Feb-2023.mat')
% load('BATS_12_iteration_10000wEMGnorm_TiedPostRamp13-Mar-2023')
% figure
% subplot(2,1,1)
% hold on
%
% early_reactive1=nanmean(X1(:,eA),2);
% late_reactive1=nanmean(X1(:,lA),2);
% reactive1=nanmean(X1(:,eP),2);
% late_Preactive1=nanmean(X1(:,lP),2);
% histogram(reactive1)
% xlabel('Reactive')
% disp('BATS - Reactive')
%
% [h,p,ci,stats] = ttest(reactive1)
%
% subplot(2,1,2)
% hold on
% context1=nanmean(X2(:,eP),2);
% early_context1=nanmean(X2(:,eA),2);
% late_context1=nanmean(X2(:,lA),2);
% late_Pcontext1=nanmean(X2(:,lP),2);
% histogram(context1)
% disp('BATS - Context')
% [h,p,ci,stats]= ttest(context1)
% xlabel('Context')
%
%
% % early_reactive2=nanmean(X3(:,eA),2);
% % late_reactive2=nanmean(X3(:,lA),2);
% % reactive2=nanmean(X3(:,eP),2);
% % late_Preactive2=nanmean(X3(:,lP),2);
%
% load('BATR_12_iteration_10000wEMGnorm_TiedPostRamp13-Mar-2023')
% % load('BATR_12_iteration_10000wEMGnorm_16-Feb-2023.mat')
% % load('BATR_11_Individual_C.mat')
% figure
% subplot(2,1,1)
% hold on
% reactive=nanmean(X1(:,eP),2);
% early_reactive=nanmean(X1(:,eA),2);
% late_reactive=nanmean(X1(:,lA),2);
% late_Preactive=nanmean(X1(:,lP),2);
% %
% % reactive3=nanmean(X3(:,eP),2);
% % early_reactive3=nanmean(X3(:,eA),2);
% % late_reactive3=nanmean(X3(:,lA),2);
% % late_Preactive3=nanmean(X3(:,lP),2);
%
% histogram(reactive)
% xlabel('Reactive')
% disp('BATR - Reactive')
% [h,p,ci,stats] = ttest(reactive)
% subplot(2,1,2)
% hold on
% context=nanmean(X2(:,eP),2);
% early_context=nanmean(X2(:,eA),2);
% late_context=nanmean(X2(:,lA),2);
% late_Pcontext=nanmean(X2(:,lP),2);
% histogram(context)
% disp('BATR - Context')
% [h,p,ci,stats]= ttest(context)
% xlabel('Context')
% legend('BATS','BATR')
% set(gcf,'color','w')
%
%
% figure
% % set(gca,'TickLabelInterpreter','latex');
%
% % stats_eP=modifiedBoxPlot([1:6],[reactive1 reactive reactive2 reactive3 context1 context]);
% % set(gca,'XTick',[1 2,3,4,5,6],'XTickLabel',{'W^{OG}_{intro}','W^{TM}_{intro}','W^{OG}_{removal}','W^{TM}_{removal}','W^{OG}_{contex}','W^{TM}_{contex}'},'FontSize',10)
% stats_eP=modifiedBoxPlot([1:4],[-reactive1 -reactive context1 context]);
% set(gca,'XTick',[1 2,3,4],'XTickLabel',{'W^{OG}_{intro}','W^{TM}_{intro}','W^{OG}_{contex}','W^{TM}_{contex}'},'FontSize',10)
%
% title('Early Post-adaptation')
% stats_eP.actualCI(1,:)=stats_eP.mean-abs(stats_eP.CI(1,:));
% stats_eP.actualCI(2,:)=stats_eP.mean+abs(stats_eP.CI(2,:));
%
% display('Reactive Distribution test')
% [h,p,ks2stat] = kstest2(reactive1,reactive)
%
%
% display('Context Distribution test')
% [h,p,ks2stat] = kstest2(context1, context)
%
%
% % stats_lP=modifiedBoxPlot([1:6],[late_Preactive1 late_Preactive late_Preactive2 late_Preactive3 late_Pcontext1 late_Pcontext]);
% % set(gca,'XTick',[1 2,3,4,5,6],'XTickLabel',{'W^{OG}_{intro}','W^{TM}_{intro}','W^{OG}_{removal}','W^{TM}_{removal}','W^{TM}_{contex}'},'FontSize',10)
% stats_lP=modifiedBoxPlot([1:4],[-late_Preactive1 -late_Preactive late_Pcontext1 late_Pcontext]);
% set(gca,'XTick',[1 2,3,4,],'XTickLabel',{'W^{OG}_{intro}','W^{TM}_{intro}','W^{TM}_{contex}'},'FontSize',10)
%
% title('Late Post-adaptation')
%
% stats_lP.actualCI(1,:)=stats_lP.mean-abs(stats_lP.CI(1,:));
% stats_lP.actualCI(2,:)=stats_lP.mean+abs(stats_lP.CI(2,:));
%
% %% Expert behavior
%
% lE=1030-45:1030-5;
%
% % load('BATS_12_iteration_10000wEMGnorm_16-Feb-2023.mat')
% figure
% subplot(2,1,1)
% hold on
%
%
% lateExpert_reactive=nanmean(X1(:,lE),2);
% histogram(lateExpert_reactive)
% xlabel('Reactive')
% disp('BAT - Reactive')
%
% subplot(2,1,2)
% hold on
%
% lateExpert_Context=nanmean(X2(:,lE),2);
% histogram(lateExpert_Context)
% disp('BAT - Context')
%
% xlabel('Context')
%
%
% set(gcf,'color','w')
%
%
%
% figure
% % set(gca,'TickLabelInterpreter','latex');
%
% stats_eP=modifiedBoxPlot([1:2],[lateExpert_reactive lateExpert_Context]);
% set(gca,'XTick',[1 2],'XTickLabel',{'LateReactive_{expert}','LateContext_{expert}'},'FontSize',10)
% title('Expert')
%
%
% %%
% stats_OG=modifiedBoxPlot([1:4],[reactive1 context1 late_Preactive1 late_Pcontext1]);
% set(gca,'XTick',[1 2,3,4],'XTickLabel',{'W^{OG}_{earlyP-reactive}','W^{OG}_{earlyP-reactive}','W^{OG}_{lateP-contex}','W^{OG}_{lateP-contex}'},'FontSize',10)
% title('OG Post-adaptation')
% ylim([-1 1.5])
%
%
% stats_TM=modifiedBoxPlot([1:4],[reactive context late_Preactive late_Pcontext]);
% set(gca,'XTick',[1 2,3,4],'XTickLabel',{'W^{TM}_{earlyP-reactive}','W^{TM}_{earlyP-reactive}','W^{TM}_{lateP-contex}','W^{T<}_{lateP-contex}'},'FontSize',10)
% title('TM Post-adaptation')
% ylim([-1 1.5])
% subplot(2,1,1)
% hold on
%
%
% lateExpert_reactive=nanmean(X1(:,lE),2);
% histogram(lateExpert_reactive)
% xlabel('Reactive')
% disp('BAT - Reactive')
%
% subplot(2,1,2)
% hold on
%
% lateExpert_Context=nanmean(X2(:,lE),2);
% histogram(lateExpert_Context)
% disp('BAT - Context')
%
% xlabel('Context')
%
%
% set(gcf,'color','w')
%
%
%
% figure
% % set(gca,'TickLabelInterpreter','latex');
%
% stats_eP=modifiedBoxPlot([1:2],[lateExpert_reactive lateExpert_Context]);
% set(gca,'XTick',[1 2],'XTickLabel',{'LateReactive_{expert}','LateContext_{expert}'},'FontSize',10)
% title('Expert')
%
%
% %%
% stats_OG=modifiedBoxPlot([1:4],[reactive1 context1 late_Preactive1 late_Pcontext1]);
% set(gca,'XTick',[1 2,3,4],'XTickLabel',{'W^{OG}_{earlyP-reactive}','W^{OG}_{earlyP-reactive}','W^{OG}_{lateP-contex}','W^{OG}_{lateP-contex}'},'FontSize',10)
% title('OG Post-adaptation')
% ylim([-1 1.5])
%
%
% stats_TM=modifiedBoxPlot([1:4],[reactive context late_Preactive late_Pcontext]);
% set(gca,'XTick',[1 2,3,4],'XTickLabel',{'W^{TM}_{earlyP-reactive}','W^{TM}_{earlyP-reactive}','W^{TM}_{lateP-contex}','W^{T<}_{lateP-contex}'},'FontSize',10)
% title('TM Post-adaptation')
% ylim([-1 1.5])
