% This is a script to get the confidance interval of the step-by-step
% weights

%%
% upload your path 
addpath(genpath('/Users/dulcemariscal/Documents/GitHub/Generalization_Regressions'))
addpath(genpath('/Users/dulcemariscal/Documents/GitHub/labTools'))
addpath(genpath('/Users/dulcemariscal/Documents/GitHub/LongAdaptation'))
addpath(genpath('/Users/dulcemariscal/Documents/GitHub/R01'))
addpath(genpath('/Users/dulcemzariscal/Documents/GitHub/splitbelt-EMG-adaptation'))
addpath(genpath('/Users/dulcemariscal/Documents/GitHub/EMG-LTI-SSM'))
addpath(genpath('/Users/dulcemariscal/Documents/GitHub/matlab-linsys'))
rmpath(genpath('/Users/dulcezmariscal/Documents/GitHub/PittSMLlab'))
 

%% Load data and Plot checkerboard for all conditions.
clear all; close all; clc;
% clear all; clc;

groupID ='BAT';
[normalizedGroupData, newLabelPrefix,n_subjects]=creatingGroupdataWnormalizedEMG(groupID);

%% Removing bad muscles
%This script make sure that we always remove the same muscle for the
%different analysis
normalizedGroupData= RemovingBadMuscleToSubj(normalizedGroupData);

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
%     load BATS_12_IndvLegsC17_ShortPertubations_RemovedBadMuscle_1RemovBias_0.mat
elseif  strcmp(groupID,'BATR')
    fname='dynamicsData_BATR_subj_12_RemoveBadMuscles1_splits_0_WithPost2V2.h5'
%     load BATR_12_IndvLegsC16_ShortPertubations_RemovedBadMuscle_1RemovBias_0.mat
end

EMGdata2=h5read(fname,'/EMGdata');
 
binwith=10;
[~,~,~,~,~,~,EMGdata,labels]=groupDataToMatrixForm_Update(1:size(EMGdata2,3),fname,0);
muscPhaseIdx=1:size(EMGdata,2); %

% context= find(strcmp(epochOfInterest,'Optimal')==1);
% reactive2=find(strcmp(epochOfInterest,'NegShort_{late}')==1);



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
%             temp3(:,:,i,l)=pinv(unit(:,:,i,l)'); %geeting the inverse of the constant
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

%% Group data plotting 
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
          Xhat(:,:,i,l) =temp3(:,:,i,l)'*Ymuscles(:,:,i)'; %x= y/C
            Yhat(:,:,i,l)=  unit(:,:,i,l)* Xhat(:,:,i,l) ; %Estimated Y with the constants 
            dynamics(:,:,i,l)=Xhat(:,:,i,l)'; %step-by-step dynamics
            
            
       
        
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

%% Group data plotting 
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
        
        x0=x_mean; % x0,y0 ellipse centre coordinates
        y0=y_mean;
 
        
        text(x0+.02,y0,{labels(m).Data(1:end-1)})
        if P_y(1)<0 &&  P_y(2)>0 %P_TM(1)<0 &&  P_TM(2)>0 ||
            rectangle('Position',[llc,CIrng],'Curvature',[1,1],'EdgeColor',clrMap(m+3,:),'LineStyle','--');
        else
            rectangle('Position',[llc,CIrng],'Curvature',[1,1],'EdgeColor',clrMap(m+3,:));
        end

        
        plot(centers(1), centers(2), 'o', 'MarkerFaceColor', clrMap(m+3,:), 'MarkerSize',10, 'LineWidth', 1,'MarkerEdgeColor',clrMap(m+3,:))%clrMap(m+3,:))
        xlabel('Treadmill')
        ylabel('Overground')

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
  
end

%% Individual muscle data plotting 
clrMap = colorcube(28*3);
muscles=[1:14 1:14];
g=[14:-1:1 14:-1:1];
ff=[1:14 1:14;15:28 15:28];
range=481:485;
for dyn=2%

    temp=[];
    x=[];
    
    for m=1:28
        figure(ff(dyn,m))
        hold on
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
        
        x0=x_mean; % x0,y0 ellipse centre coordinates
        y0=y_mean;

        
        text(x0+.02,y0,{labels(m).Data(1:end-1)})
        if P_y(1)<0 &&  P_y(2)>0 || P_x(1)<0 &&  P_x(2)>0 
            rectangle('Position',[llc,CIrng],'Curvature',[1,1],'EdgeColor',clrMap(m+3,:),'LineStyle','--');
        else
            rectangle('Position',[llc,CIrng],'Curvature',[1,1],'EdgeColor',clrMap(m+3,:));
        end

        
        plot(centers(1), centers(2), 'o', 'MarkerFaceColor', clrMap(m+3,:), 'MarkerSize',10, 'LineWidth', 1,'MarkerEdgeColor',clrMap(m+3,:))%clrMap(m+3,:))
        xlabel('Treadmill')
        ylabel('Overground')
%         
        if m<15
            Li{1}=scatter(nanmean(TM_2(dyn,range,m)),nanmean(OG_2(dyn,range,m)),100,"filled",'MarkerFaceColor', 'b');
            text(nanmean(TM_2(dyn,range,m))+.02,nanmean(OG_2(dyn,range,m)),{labels(m).Data(1:end-1)})
        else
            Li{2}=scatter(nanmean(TM_2(dyn,range,m)),nanmean(OG_2(dyn,range,m)),100,"filled",'MarkerFaceColor', 'r')  ;
               text(nanmean(TM_2(dyn,range,m))+.02,nanmean(OG_2(dyn,range,m)),{labels(m).Data(1:end-1)})
        end

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
    xlim([-.5 1])
    ylim([-.5 1])
end
xlabel('Treadmill')
ylabel('Overground')
    set(gcf,'color','w')
    end

end
