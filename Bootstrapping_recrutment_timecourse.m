% This is a script to get the confidance interval of the step-by-step
% weights 

%% Load data and Plot checkerboard for all conditions.
% clear all; close all; clc;
clear all; clc;

groupID ='BATS';
[normalizedGroupData, newLabelPrefix,n_subjects]=creatingGroupdataWnormalizedEMG(groupID);
%% Removing bad muscles 
%This script make sure that we always remove the same muscle for the
%different analysis 
normalizedGroupData= RemovingBadMuscleToSubj(normalizedGroupData);
%%  Getting the step-by-step data 

%Adaptation epochs
strides=[-40 440 200];
cond={'TM base','Adaptation','Post 1'}; %Conditions for this group

exemptFirst=[1];
exemptLast=[5]; %Strides needed
names={};
shortNames={};

ep=defineEpochs(cond,cond,strides,exemptFirst,exemptLast,'nanmean',{'Base','Adapt','Post1'}); %epochs

padWithNaNFlag=true;
[dataEMG,labels,allDataEMG]=normalizedGroupData.getPrefixedEpochData(newLabelPrefix(end:-1:1),ep,padWithNaNFlag); %Getting the data 

%Flipping EMG:
for i=1:length(allDataEMG)
    aux=reshape(allDataEMG{i},size(allDataEMG{i},1),size(labels,1),size(labels,2),size(allDataEMG{i},3));
    allDataEMG{i}=reshape(flipEMGdata(aux,2,3),size(aux,1),numel(labels),size(aux,4));
end

EMGdata=cell2mat(allDataEMG);
muscPhaseIdx=1:size(EMGdata,2);


%% Getting the C values 
epochOfInterest={'TM base','TM mid 1','PosShort_{early}','PosShort_{late}','Ramp','Optimal','Adaptation','Adaptation_{early}','TiedPostPos','TMmid2','NegShort_{late}','Post1_{Early}','TMbase_{early}'};
ep=defineRegressorsDynamicsFeedback('nanmean');
refEpTM = defineReferenceEpoch('TM base',ep);
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

%% Bootstrapping 

bootstrap=1 %Do you want to run the loop (1 yes 0 No) 
X1=[];
X2=[];
replacement=1 %do you want to do it with replacement (1 yes 0 No) 


if bootstrap
    if replacement
        n=1000; %number of iterations
    else
        n=1
    end
    
    f = waitbar(0,'1','Name','Boostrapping Data',...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    
    setappdata(f,'canceling',0);
    
    for i=1:n
        
        temp=[];
        temp2=[];
        DataBoot={};
        % Check for clicked Cancel button
        if getappdata(f,'canceling')
            break
        end
        % Update waitbar and message
        ww=waitbar(i/n,f,['Iteration ' num2str(i)]);
        
        %         normalizedGroupData =[];
        if replacement
            subjIdx=datasample(1:n_subjects,n_subjects,'Replace',true);
        else
            subjIdx=datasample(1:n_subjects,n_subjects,'Replace',false);
        end
        DataBoot=Data(subjIdx,:);
        
        for c=1:length(Data)
            for s=1:n_subjects
                temp(:,:,s)=DataBoot{s,c};
            end
            temp2{1,c}=nanmedian(temp,3); 
            tt=temp2{1,c}(:,end:-1:1);
            temp2{1,c}=tt;
            temp2{1,c}=reshape(temp2{1,c},14*2*12,1);
        end
        
        temp2=cell2mat(temp2)';
        temp2=temp2-fftshift(temp2,2);
        temp2=temp2(:,1:size(temp2,2)/2,:);
        temp2=temp2';
        C2=[temp2(:,5) temp2(:,6)];
        
        Y2=EMGdata(:,muscPhaseIdx,subjIdx);
        Y2=nanmedian(Y2,3); 
        Yasym2=Y2-fftshift(Y2,2);
        Yasym2=Yasym2(:,1:size(Yasym2,2)/2,:);
        
        
        
        
        Cinv2=pinv(C2);
        Xdynamics= Cinv2*Yasym2'; %x= y/C
     
        X1=[X1;Xdynamics(1,:)];
        X2=[X2;Xdynamics(2,:)];
        

        
        
    end
    close(ww)
    
end

delete(f)

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

d{1} = X1(:,1:40);
d{2} = X1(:,41:480);
d{3} =  X1(:,481:end);
index{1}= 1:40;
index{2}= 42:481;
index{3}= 483:682;

for i=1:3
    y = nanmean(d{i},1); % your mean vector;
    %     x = 1:numel(y);
    x=index{i};
    std_dev = nanstd(d{i},1);
    curve1 = y + std_dev;
    curve2 = y - std_dev;
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    fill(x2, inBetween,'r','FaceAlpha',0.3,'EdgeColor','none')
    hold on;
    ylabel('W_{reactive}')
    xlabel('Strides')
    plot(x, y, 'r', 'LineWidth', 2);
   
end
yline(0)
% legend('Training')


subplot(2,1,2)
hold on
% X_2=[X2(:,1:40) nan(size(X2(:,1:40),1),1) X2(:,41:480) nan(size(X2(:,1:40),1),1) X2(:,481:end)];
% y = nanmean(X2,1); % your mean vector;
% x = 1:numel(y);
% std_dev = nanstd(X2,1);
% curve1 = y + std_dev;
% curve2 = y - std_dev;
% x2 = [x, fliplr(x)];
% inBetween = [curve1, fliplr(curve2)];
% fill(x2, inBetween,'r','FaceAlpha',.3,'EdgeColor','none')
% hold on;
% ylabel('W_{context}')
% xlabel('Strides')
% plot(x, y, 'r', 'LineWidth', 2);
d{1} = X2(:,1:40);
d{2} = X2(:,41:480);
d{3} =  X2(:,481:end);

for i=1:3
    y = nanmean(d{i},1); % your mean vector;
    %     x = 1:numel(y);
    x=index{i};
    std_dev = nanstd(d{i},1);
    curve1 = y + std_dev;
    curve2 = y - std_dev;
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    fill(x2, inBetween,'r','FaceAlpha',0.3,'EdgeColor','none')
    hold on;
    ylabel('W_{context}')
    xlabel('Strides')
    plot(x, y, 'r', 'LineWidth', 2);
    
end
yline(0)
set(gcf,'color','w')
%%
save([groupID,'_',num2str(n_subjects),'_iteration_', num2str(n)],'X1','X2','subID')


%% EMG norm 
load('BATR_12_iteration_10000wEMGnorm_16-Feb-2023.mat')
subplot(2,1,1)
hold on
TR_post=nanmean(EMGnorm(:,481:485),2);
TR_eA=nanmean(EMGnorm(:,41:51),2);
TR_lA=nanmean(EMGnorm(:,480-40:480),2);
histogram(TR_post,'FaceColor','r')
xlabel('EMGnorm')
disp('BATR')
% [h,p,ci,stats] = ttest(reactive)
% subplot(2,1,2)
hold on

load('BATS_12_iteration_10000wEMGnorm_16-Feb-2023.mat')
TS_post=nanmean(EMGnorm(:,481:485),2);
TS_eA=nanmean(EMGnorm(:,41:51),2);
TS_lA=nanmean(EMGnorm(:,480-40:480),2);
histogram(TS_post,'FaceColor','b')
disp('BATR - Context')
% [h,p,ci,stats]= ttest(context)
xlabel('Context')
legend('BATR','BATS')
set(gcf,'color','w')

disp('Context Distribution test')
[h,p,ks2stat] = kstest2( TS_post, TR_post)

% disp('t-test')
% [h,p,ks2stat] = ttest2( TS_post, TR_post)

title(['Two-sample Kolmogorov-Smirnov, p=',num2str(p),' Rejects the null hypothesis =', num2str(h)])

stats_norm=modifiedBoxPlot([1,2],[ TS_post TR_post]);
% set(gca,'title','Marce')
set(gca,'XTick',[1 2],'XTickLabel',{'EMGnorm^{TS}_{post}','EMGnorm^{TR}_{post}' },'FontSize',10)
title('Post-adaptation')

stats_norm=modifiedBoxPlot([1,2,3,4],[ TS_eA TR_eA TS_lA TR_lA]);
% set(gca,'title','Marce')
set(gca,'XTick',[1 2,3,4],'XTickLabel',{'EMGnorm^{TS}_{early}','EMGnorm^{TR}_{early}' ,'EMGnorm^{TS}_{late}','EMGnorm^{TR}_{late}' },'FontSize',10)
title('Adaptation')

%% Adaptation 

eA=42:51;
lA=435:480-5;
eP=482:486;


load('BAT_24_iteration_10000wEMGnorm_16-Feb-2023.mat')
figure
subplot(2,1,1)
hold on
early_reactive1=nanmean(X1(:,eA),2);
late_reactive1=nanmean(X1(:,lA),2);


early_context1=nanmean(X2(:,eA),2);
late_context1=nanmean(X2(:,lA),2);
histogram(early_context1)


stats_Adapt=modifiedBoxPlot([1:4],[early_reactive1  early_context1   late_reactive1 late_context1 ]);
set(gca,'XTick',[1 2,3,4],'XTickLabel',{'W_{early-reactive}','W_{early-contex}','W_{late-reactive}','W_{late-contex}'},'FontSize',10)
title('Adaptation (All Data)')

%% Post-adaptation data 
eA=42:51;
lA=435:480-5;
eP=482:486;
lP=680-45:680-5;

% load('BATS_12_iteration_10000wEMGnorm_16-Feb-2023.mat')
load('BATS_12_iteration_10000wEMGnorm_TiedPostRamp13-Mar-2023')
figure
subplot(2,1,1)
hold on

early_reactive1=nanmean(X1(:,eA),2);
late_reactive1=nanmean(X1(:,lA),2);
reactive1=nanmean(X1(:,eP),2);
late_Preactive1=nanmean(X1(:,lP),2);
histogram(reactive1)
xlabel('Reactive')
disp('BATS - Reactive')

[h,p,ci,stats] = ttest(reactive1)

subplot(2,1,2)
hold on
context1=nanmean(X2(:,eP),2);
early_context1=nanmean(X2(:,eA),2);
late_context1=nanmean(X2(:,lA),2);
late_Pcontext1=nanmean(X2(:,lP),2);
histogram(context1)
disp('BATS - Context')
[h,p,ci,stats]= ttest(context1)
xlabel('Context')


% early_reactive2=nanmean(X3(:,eA),2);
% late_reactive2=nanmean(X3(:,lA),2);
% reactive2=nanmean(X3(:,eP),2);
% late_Preactive2=nanmean(X3(:,lP),2);

load('BATR_12_iteration_10000wEMGnorm_TiedPostRamp13-Mar-2023')
% load('BATR_12_iteration_10000wEMGnorm_16-Feb-2023.mat')
% load('BATR_11_Individual_C.mat')
figure
subplot(2,1,1)
hold on
reactive=nanmean(X1(:,eP),2);
early_reactive=nanmean(X1(:,eA),2);
late_reactive=nanmean(X1(:,lA),2);
late_Preactive=nanmean(X1(:,lP),2);
% 
% reactive3=nanmean(X3(:,eP),2);
% early_reactive3=nanmean(X3(:,eA),2);
% late_reactive3=nanmean(X3(:,lA),2);
% late_Preactive3=nanmean(X3(:,lP),2);

histogram(reactive)
xlabel('Reactive')
disp('BATR - Reactive')
[h,p,ci,stats] = ttest(reactive)
subplot(2,1,2)
hold on
context=nanmean(X2(:,eP),2);
early_context=nanmean(X2(:,eA),2);
late_context=nanmean(X2(:,lA),2);
late_Pcontext=nanmean(X2(:,lP),2);
histogram(context)
disp('BATR - Context')
[h,p,ci,stats]= ttest(context)
xlabel('Context')
legend('BATS','BATR')
set(gcf,'color','w')


figure 
% set(gca,'TickLabelInterpreter','latex');

% stats_eP=modifiedBoxPlot([1:6],[reactive1 reactive reactive2 reactive3 context1 context]);
% set(gca,'XTick',[1 2,3,4,5,6],'XTickLabel',{'W^{OG}_{intro}','W^{TM}_{intro}','W^{OG}_{removal}','W^{TM}_{removal}','W^{OG}_{contex}','W^{TM}_{contex}'},'FontSize',10)
stats_eP=modifiedBoxPlot([1:4],[-reactive1 -reactive context1 context]);
set(gca,'XTick',[1 2,3,4],'XTickLabel',{'W^{OG}_{intro}','W^{TM}_{intro}','W^{OG}_{contex}','W^{TM}_{contex}'},'FontSize',10)

title('Early Post-adaptation')
stats_eP.actualCI(1,:)=stats_eP.mean-abs(stats_eP.CI(1,:));
stats_eP.actualCI(2,:)=stats_eP.mean+abs(stats_eP.CI(2,:));

display('Reactive Distribution test')
[h,p,ks2stat] = kstest2(reactive1,reactive)


display('Context Distribution test')
[h,p,ks2stat] = kstest2(context1, context)


% stats_lP=modifiedBoxPlot([1:6],[late_Preactive1 late_Preactive late_Preactive2 late_Preactive3 late_Pcontext1 late_Pcontext]);
% set(gca,'XTick',[1 2,3,4,5,6],'XTickLabel',{'W^{OG}_{intro}','W^{TM}_{intro}','W^{OG}_{removal}','W^{TM}_{removal}','W^{TM}_{contex}'},'FontSize',10)
stats_lP=modifiedBoxPlot([1:4],[-late_Preactive1 -late_Preactive late_Pcontext1 late_Pcontext]);
set(gca,'XTick',[1 2,3,4,],'XTickLabel',{'W^{OG}_{intro}','W^{TM}_{intro}','W^{TM}_{contex}'},'FontSize',10)

title('Late Post-adaptation')

stats_lP.actualCI(1,:)=stats_lP.mean-abs(stats_lP.CI(1,:));
stats_lP.actualCI(2,:)=stats_lP.mean+abs(stats_lP.CI(2,:));

%% Expert behavior 

lE=1030-45:1030-5;

% load('BATS_12_iteration_10000wEMGnorm_16-Feb-2023.mat')
figure
subplot(2,1,1)
hold on


lateExpert_reactive=nanmean(X1(:,lE),2);
histogram(lateExpert_reactive)
xlabel('Reactive')
disp('BAT - Reactive')

subplot(2,1,2)
hold on

lateExpert_Context=nanmean(X2(:,lE),2);
histogram(lateExpert_Context)
disp('BAT - Context')

xlabel('Context')


set(gcf,'color','w')



figure 
% set(gca,'TickLabelInterpreter','latex');

stats_eP=modifiedBoxPlot([1:2],[lateExpert_reactive lateExpert_Context]);
set(gca,'XTick',[1 2],'XTickLabel',{'LateReactive_{expert}','LateContext_{expert}'},'FontSize',10)
title('Expert')


%%
stats_OG=modifiedBoxPlot([1:4],[reactive1 context1 late_Preactive1 late_Pcontext1]);
set(gca,'XTick',[1 2,3,4],'XTickLabel',{'W^{OG}_{earlyP-reactive}','W^{OG}_{earlyP-reactive}','W^{OG}_{lateP-contex}','W^{OG}_{lateP-contex}'},'FontSize',10)
title('OG Post-adaptation')
ylim([-1 1.5])


stats_TM=modifiedBoxPlot([1:4],[reactive context late_Preactive late_Pcontext]);
set(gca,'XTick',[1 2,3,4],'XTickLabel',{'W^{TM}_{earlyP-reactive}','W^{TM}_{earlyP-reactive}','W^{TM}_{lateP-contex}','W^{T<}_{lateP-contex}'},'FontSize',10)
title('TM Post-adaptation')
ylim([-1 1.5])
