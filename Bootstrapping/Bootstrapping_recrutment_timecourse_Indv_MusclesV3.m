% This is a script to get the confidance interval of the step-by-step
% weights. The regressors are the same for each iteration and group

%%
% upload your path
main='/Users/dulcemariscal/Documents/GitHub/';
% main='C:\Users\dum5\OneDrive - University of Pittsburgh\_BoxMigration\GitHub\';
addpath(genpath([main, 'Generalization_Regressions']))
addpath(genpath([main,'labTools']))
addpath(genpath([main,'R01']))
addpath(genpath([main,'splitbelt-EMG-adaptation']))
addpath(genpath([main,'EMG-LTI-SSM']))
addpath(genpath([main,'matlab-linsys']))
% rmpath(genpath([main,'PittSMLlab']))

%% Load data and Plot checkerboard for all conditions.
clear all; close all; clc;
% clear all; clc;

adaptationonly=0; 
removeMuscles=1;

if adaptationonly==1
    groupID ={'BAT'};
else
    groupID ={'BATR','BATS'};
end

EMG2=[];
Regres=[];
r2=[];
r2uncenter=[];

  
rng('default') % for reproducibily we set the random seed 

for id=1:length(groupID)

    [normalizedGroupData, newLabelPrefix,n_subjects]=creatingGroupdataWnormalizedEMG(groupID{id});

    %% Removing bad muscles
    %This script make sure that we always remove the same muscle for the
    %different analysis
    if removeMuscles==1
        normalizedGroupData= RemovingBadMuscleToSubj(normalizedGroupData);
    end
%% Getting the regressors values
%% 
ep=defineRegressorsDynamicsFeedback('nanmean');
epochOfInterest={'TM base','NegShort_{late}','Ramp','Optimal'};

padWithNaNFlag=true; %If no enough steps fill with nan, let this on

for l=1:length(epochOfInterest)
    
    ep2=defineReferenceEpoch(epochOfInterest{l},ep);
    
    [dataEMG,labels,allDataEMG2]=normalizedGroupData.getPrefixedEpochData(newLabelPrefix(end:-1:1),ep2,padWithNaNFlag); %Getting the data
    
    %Flipping EMG:0
    for i=1:length(allDataEMG2)
        aux=reshape(allDataEMG2{i},size(allDataEMG2{i},1),size(labels,1),size(labels,2),size(allDataEMG2{i},3));
        allDataEMG2{i}=reshape(flipEMGdata(aux,2,3),size(aux,1),numel(labels),size(aux,4));
        
    end
    
    Data{l}=squeeze(nanmean(allDataEMG2{:},1));    
end
    Regres=[Regres ; Data];
    %%  Getting the step-by-step data
    % Adaptation epochs

    if adaptationonly==1
        strides=[-40 440];
        cond={'TM base','Adaptation'}; %Conditions for this group
    else
        strides=[-40 200];
        if contains(groupID{id},'TR') %for treadmill Post 1
            cond={'TM base','Post 1'}; %Conditions for this group
        else % for overground post 1
            cond={'OG base','Post 1'}; %Conditions for this group
        end
    end

    exemptFirst=[1];
    exemptLast=[5]; %Strides needed
    names={};
    shortNames={};
    if adaptationonly==1
        ep=defineEpochs(cond,cond,strides,exemptFirst,exemptLast,'nanmean',{'Base','Adapt'}); %epochs
    else
        ep=defineEpochs(cond,cond,strides,exemptFirst,exemptLast,'nanmean',{'Base','Post1'}); %epochs

    end
    padWithNaNFlag=true; %If no enough steps fill with nan, let this on
    [dataEMG,labels,allDataEMG]=normalizedGroupData.getPrefixedEpochData(newLabelPrefix(end:-1:1),ep,padWithNaNFlag); %Getting the data

    %Flipping EMG:0
    for i=1:length(allDataEMG)
        aux=reshape(allDataEMG{i},size(allDataEMG{i},1),size(labels,1),size(labels,2),size(allDataEMG{i},3));
        allDataEMG{i}=reshape(flipEMGdata(aux,2,3),size(aux,1),numel(labels),size(aux,4));
    end

    EMGdata=cell2mat(allDataEMG); %Getting EMG data per participants

    EMG2=cat(3, EMG2, EMGdata);
end
%% Bootstrapping
%% Unnecesary but idk how to fix it 
if adaptationonly==0

    for  l=1:length(epochOfInterest)
        Regres{l}=cell2mat([Regres(1,l),Regres(2,l)])
    end
    
end

%%
muscPhaseIdx=1:size(EMGdata,2); %
epochOfInterest={'TM base','NegShort_{late}','Ramp','Optimal'};
context= find(strcmp(epochOfInterest,'Optimal')==1);

if adaptationonly==1
    reactive=find(strcmp(epochOfInterest,'Ramp')==1);
else
    reactive=find(strcmp(epochOfInterest,'NegShort_{late}')==1);
end

base=find(strcmp(epochOfInterest,'TM base')==1);
% 
matrix=[];

%Data TR
if length(groupID)>1
    TR=EMG2(:,:,1:12);
    TS=EMG2(:,:,13:24);
else
    TR=EMG2;
end

bootstrap=1; %Do you want to run the loop (1 yes 0 No)
X1=[];
X2=[];
replacement=0; %do you want to do it with replacement (1 yes 0 No)
regre_Const=1; % To keep the regressors constants for both groups

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
    C_inv=nan(12,2,28,n);
    Yhat=nan(12,size(EMG2,1),28,n);
    temp4=nan(size(EMG2,1),2,28,n);

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

        if  regre_Const
            subjIdx=datasample(1:24,24,'Replace',false);
        else
            subjIdx=datasample(1:24,24,'Replace',true);
        end

        if replacement %doing the bootstrap with replacement
            if adaptationonly==1
                groupIdx=datasample(1:24,24,'Replace',true);
            else
                groupIdx=datasample(1:12,12,'Replace',true);
            end
        else
            groupIdx=datasample(1:12,12,'Replace',false);
        end

%         DataBoot=Data(:,subjIdx); %Subject pick at each loop fro regressors

        %This loop is to compute the constrant on our regressions
        if l==1 && regre_Const==1
            for c=1:length(epochOfInterest)
                    temp= Regres{c}(:,subjIdx); %DataBoot{s,c};

                x{1,c}=nanmedian(temp,2);
          end


            x=cell2mat(x)';
            x=x';

            %C values that we are using for the regressions
%             C2=[x(:,reactive) x(:,context)]- x(:,base);
            C2=[x(:,reactive) x(:,context)]- x(:,base);

            %reorganize the data to be separatend by muscle
            Cmuscles=reshape(C2',2,12,28);

        end
        %Picking the data muscles that we want and participants
        Y_TR=TR(:,muscPhaseIdx,groupIdx);

        if length(groupID)>1
            Y_TS=TS(:,muscPhaseIdx,groupIdx);
        end

        %removing the bias for group
        Y_TR=nanmedian(Y_TR,3); %getting the median of the group
        bias=nanmean(Y_TR(1:40,:,:)) ; %estimating the gorup baseline
        Y_TR=Y_TR-bias; %removing the bias from the data


        if length(groupID)>1
            Y_TS=nanmedian(Y_TS,3); %getting the median of the group
            bias=nanmean(Y_TS(1:40,:,:)) ; %estimating the group baseline
            Y_TS=Y_TS-bias; %removing the bias from the data
        end

        %reorganize the data to be separatend by muscle
        Ymuscles_TR=reshape(Y_TR(1:size(EMG2,1),:),size(EMG2,1),12,28);

        if length(groupID)>1
            Ymuscles_TS=reshape(Y_TS(1:size(EMG2,1),:),size(EMG2,1),12,28);
        end

        %Linear regression individual muscles
        reconstruction_indv=[];
        data=[];
        C_indv=[];
        X_indv=[];
        data=[];
        C_indv=[];
        reconstruction_indv=[];
        reconstruction_contextual=[];
        reconstruction_reactive=[];
        data_shifted =[];

        for i=1:size(Ymuscles_TR,3) %loop for individual muscle fit

            % Getting the inverse of the regressors
            unit(:,:,i)=Cmuscles(:,:,i)'./vecnorm(Cmuscles(:,:,i)');
            if l==1            
                matrix=[matrix;unit(:,:,i)];
            end
             C_inv=pinv(unit(:,:,i)'); %geting the inverse of the constant

            %%% TR
            Xhat_TR(:,:,i,l) =C_inv'*Ymuscles_TR(:,:,i)'; %x= y/C
            Yhat_TR(:,:,i,l)=  unit(:,:,i)* Xhat_TR(:,:,i,l) ; %Estimated Y with the constants
            dynamics_TR(:,:,i,l)=Xhat_TR(:,:,i,l)'; %step-by-step dynamics

            if length(groupID)>1
                %%% TS
                Xhat_TS(:,:,i,l) =C_inv'*Ymuscles_TS(:,:,i)'; %x= y/C
                Yhat_TS(:,:,i,l)=  unit(:,:,i)* Xhat_TS(:,:,i,l) ; %Estimated Y with the constants
                dynamics_TS(:,:,i,l)=Xhat_TS(:,:,i,l)'; %step-by-step dynamics
                
                
                
                yContextCurr = unit(:,2,i)*  Xhat_TS(2,:,i);
                yReactiveCurr = unit(:,1,i)*  Xhat_TS(1,:,i);
                
                r = abs(nanmin(yReactiveCurr)); %reactive min
                c = abs(nanmin(yContextCurr)); %contextual min
                
                reconstruction_indv =[ reconstruction_indv ; Yhat_TS(:,:,i)+c+r]; % Concatenating the data reconstructed and adding the minimium of the reactive and the contextual
                reconstruction_contextual=[reconstruction_contextual; yContextCurr+c];
                reconstruction_reactive=[reconstruction_reactive; yReactiveCurr+r];
                data_shifted =[ data_shifted ; Ymuscles_TS(:,:,i)'+c+r];
                data =[ data ; Ymuscles_TS(:,:,i)'];  % Concatenating the data
                
            end

            
        end
        
        if length(groupID)>1
            data_shifted=nanmean(data_shifted(:,41:45),2);
            data_shifted = data_shifted(~isnan(data_shifted));
            
            reconstruction_indv=nanmean(reconstruction_indv(:,41:45),2);
            reconstruction_indv=reconstruction_indv(~isnan(reconstruction_indv));
            
            reconstruction_contextual=nanmean(reconstruction_contextual(:,41:45),2);
            reconstruction_contextual=reconstruction_contextual(~isnan(reconstruction_contextual));
            
            reconstruction_reactive=nanmean(reconstruction_reactive(:,41:45),2);
            reconstruction_reactive=reconstruction_reactive(~isnan(reconstruction_reactive));
            
            r2.total(l,1)=my_Rsquared_coeff(data_shifted,reconstruction_indv,1);
            r2.contextual(l,1)=my_Rsquared_coeff(data_shifted,reconstruction_contextual,1);
            r2.reactive(l,1)=my_Rsquared_coeff(data_shifted,reconstruction_reactive,1);
            
            r2uncenter.total(l,1)=my_Rsquared_coeff(data_shifted,reconstruction_indv,0);
            r2uncenter.contextual(l,1)=my_Rsquared_coeff(data_shifted,reconstruction_contextual,0);
            r2uncenter.reactive(l,1)=my_Rsquared_coeff(data_shifted,reconstruction_reactive,0);
        end
    end
    
    close(ww)

end

delete(f)
%% Plot regressors checkerboards

muscleOrder={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT', 'VM', 'VL', 'RF', 'HIP','TFL', 'GLU'};
ytl=([strcat('f',muscleOrder) strcat('s',muscleOrder)]);  %List of muscle
ytl(end:-1:1) = ytl(:);
yt=1:length(ytl);
fs=14;
ph=figure
subplot(1,2,1)   
imagesc((reshape(matrix(:,1),12,28)'))
set(gca,'XTick',[],'YTick',yt,'YTickLabel',ytl,'FontSize',fs)
set(gca,'CLim',[-1 1]*1,'FontSize',fs); %making sure that the first plot color scheme goes [-1 1] and making the name of the labels larger


subplot(1,2,2)
imagesc((reshape(matrix(:,2),12,28)'))

color4checkerboards


% Making the figure a bit prettier
fs=14; %font size
colormap(flipud(map)) %changing the color map to the one tha defined about, we flipup the matrix bc the code does L-R and we want R-L
% c=flipud('gray');
% colormap(flipud(gray))
set(gcf,'color','w'); %setting the background white
set(gca,'XTick',[],'YTick',yt,'YTickLabel',ytl,'FontSize',fs)
colorbar %Showing the colormap bar



%%
if adaptationonly==1 && regre_Const==1

    save([groupID{1}(1:3),'_',num2str(24),'_iteration_', num2str(l),'_Individual_muscles_C_constant_Adaptation','_',datestr(now,'dd-mmmm-yyyy')],'dynamics_TR','-v7.3')

elseif regre_Const==1

    save([groupID{1}(1:3),'_',num2str(24),'_iteration_', num2str(l),'_Individual_muscles_C_constant_PostAdapt_per_group','_',datestr(now,'dd-mmmm-yyyy')],'dynamics_TR','dynamics_TS','-v7.3')

end
% save([groupID{1}(1:3),'_',num2str(24),'_iteration_', num2str(l),'_Individual_muscles_C_constant_Adaptation'],'dynamics_TR','-v7.3')
% save([groupID,'_',num2str(n_subjects),'_iteration_', num2str(n),'_Individual_muscles'],'dynamics','Yhat','Ymuscles','groupID','-v7.3')
%%
load musclesLabels.mat 
% load 'BAT_24_iteration_2000_Individual_muscles.mat'
% load BAT_24_iteration_2000_Individual_muscles_C_constant_PostAdapt_per_group.mat
% load('BATS_12_iteration_2000_Individual_muscles.mat')
% load('BAT_24_iteration_2000_Individual_muscles_C_constant_PostAdaptation_per_group_25-July-2023 08:24.mat')

OG=dynamics_TS;
% load('BATR_12_iteration_2000_Individual_muscles.mat')
TM=dynamics_TR;


% load('musclesLabels.mat')
% load('BATR_indv_muscles.mat')
% TM_2=X2asym;
% load('BATS_indv_muscles.mat')
% OG_2=X2asym;
%%
load('NCM2023_Treadmill.mat')
TM_2=X2asym;
TM_2(1,:,:)=-TM_2(1,:,:);
load('NCM2023_OG.mat')
OG_2=X2asym;
OG_2(1,:,:)=-OG_2(1,:,:);
%% Individual muscle x-axis 
c=colormap(jet(28));

TR= squeeze(mean(dynamics_TR(481:485,:,:,:),1));
TS= squeeze(mean(dynamics_TS(481:485,:,:,:),1));
TR_Reactive=squeeze(TR(1,:,:));
TR_contextual=squeeze(TR(2,:,:)); 
TS_Reactive=squeeze(TS(1,:,:));
TS_contextual=squeeze(TS(2,:,:));

%%
muscleOrder={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT',...
    'VM', 'VL', 'RF', 'HIP','TFL', 'GLU'};
ytl=([strcat('f',muscleOrder) strcat('s',muscleOrder)]);  %List of muscle
ytl(end:-1:1) = ytl(:);
yt=1:length(ytl);
post1Index = 41:45;

titles = {'Reactive','Contextual'};

for cIdx=1:2 %reactive, then context, 

    f = figure('units','normalized','outerposition',[0 0 1 0.5]);%('Position', get(0, 'Screensize'));

    for mIdx = 1:size(dynamics_TR,3) %for each muscle  muscle orders start from sGLU
       postData = nanmean(squeeze(dynamics_TR(post1Index,cIdx,mIdx,:)));
       PlotHelper.plotCI(mIdx, postData,'b',ytl{mIdx},false, true) %
    end

    title([titles{cIdx} ' (Mean, 95% from bootstrap)'])
    xlim([0 size(dynamics_TR,3)+1])
    xticks(1:size(dynamics_TR,3))
    xticklabels(ytl)
    xtickangle(45)
    PlotHelper.tightMargin(gca)

end


%% Group data plotting
clrMap = colorcube(28*3);
muscles=[1:14 1:14];
g=[14:-1:1 14:-1:1];
ff=[1:14 1:14;15:28 15:28];
range=481:485;


for i=1:28
    temp{i,1}=labels(i).Data(1:end-1);
end

label2={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT',...
    'VM', 'VL', 'RF', 'HIP','TFL', 'GLU'};

%contextual
% label2={'fLG','fMG'}; %Yes TM and Yes OG
% label2={'sTFL','sHIP','sVM','sPER','sSOL','sSEMT','sSEMB','sLG','sMG'}; %Yes TM; No OG
% label2={'fGLU','fTFL','fHIP','fRF','sRF','fVL','sVL','fVM','fSEMT','fSEMB','sTA','sBF','fBF'}; %No TM and NO OG
% label2={'sGLU','fPER','fSOL','fTA'}; %no TM and Yes OG


%reactive
% label2={'fGLU','sTA','fTA','sPER','fPER','sSOL','fSOL','sLG','fLG','fMG','sMG','fBF','sBF',...
% 	'sSEMB','sSEMT','fVM','fVL','fRF','fHIP','sHIP','fTFL'}; %Yes TM and Yes OG
% label2={'sGLU','fSEMB','fSEMT'};  %Yes TM; No OG
% label2={'sTFL'};  %No TM and NO OG
% label2={'sVM','sVL','sRF'}; %NO TM and YES OG
t=[];
t(1,:)=find(contains(temp,label2));
%%
for dyn=1:2%
    figure()
    hold on
    temp=[];
    x=[];

    for m=t
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

        if P_y(1)<0 &&  P_y(2)>0 && P_x(1)<0 &&  P_x(2)>0
            rectangle('Position',[llc,CIrng],'Curvature',[1,1],'EdgeColor',clrMap(m+3,:),'LineStyle',"-",'LineWidth',1);
        elseif P_y(1)<0 &&  P_y(2)>0
            rectangle('Position',[llc,CIrng],'Curvature',[1,1],'EdgeColor',clrMap(m+3,:),'LineStyle',"-",'LineWidth',1);
        elseif P_x(1)<0 &&  P_x(2)>0
            rectangle('Position',[llc,CIrng],'Curvature',[1,1],'EdgeColor',clrMap(m+3,:),'LineStyle','-','LineWidth',1);
        else
            rectangle('Position',[llc,CIrng],'Curvature',[1,1],'EdgeColor',clrMap(m+3,:),'LineWidth',1);
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
        %         xlim([-.5 2.5])
        %         ylim([-.5 2.5])
        yline(0)
        xline(0)
    else
        title('Contextual')
        yline(0)
        xline(0)
        xlim([-.6 1.6])
        ylim([-.6 1.6])
    end

    set(gcf,'color','w')

end

%% Individual muscle data plotting
clrMap = colorcube(28*3);
muscles=[1:14 1:14];
g=[14:-1:1 14:-1:1];
ff=[1:14 1:14;15:28 15:28];
range=440:480; %Post-adapt 481:485 Early-adapt 41:45 Late adapt 440:480
% range=

for dyn=1
    figure
    temp=[];
    x=[];

    for m=1:28
%         figure(ff(dyn,m))
        subplot(2,7,ff(1,m))
        hold on
        x=squeeze(TM(:,1,m,:));
        y=squeeze(TM(:,2,m,:));

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


        text(x0+.05,y0,{labels(m).Data(1:end-1)})

        %         if P_y(1)<0 &&  P_y(2)>0 || P_x(1)<0 &&  P_x(2)>0
        if P_y(1)<0 &&  P_y(2)>0 && P_x(1)<0 &&  P_x(2)>0
            rectangle('Position',[llc,CIrng],'Curvature',[1,1],'EdgeColor',clrMap(m+3,:),'LineStyle',":",'LineWidth',1);
        elseif P_y(1)<0 &&  P_y(2)>0
            rectangle('Position',[llc,CIrng],'Curvature',[1,1],'EdgeColor',clrMap(m+3,:),'LineStyle',"-.",'LineWidth',1);
        elseif P_x(1)<0 &&  P_x(2)>0
            rectangle('Position',[llc,CIrng],'Curvature',[1,1],'EdgeColor',clrMap(m+3,:),'LineStyle','--','LineWidth',1);
        else
            rectangle('Position',[llc,CIrng],'Curvature',[1,1],'EdgeColor',clrMap(m+3,:),'LineWidth',2);
        end


        plot(centers(1), centers(2), 'o', 'MarkerFaceColor', clrMap(m+3,:), 'MarkerSize',10, 'LineWidth', 1,'MarkerEdgeColor',clrMap(m+3,:))%clrMap(m+3,:))
%         xlabel('Treadmill')
%         ylabel('Overground')
        %
%                 if m<15
%                     Li{1}=scatter(nanmean(TM_2(dyn,range,m)),nanmean(OG_2(dyn,range,m)),100,"filled",'MarkerFaceColor', 'b');
%                     text(nanmean(TM_2(dyn,range,m))+.02,nanmean(OG_2(dyn,range,m)),{labels(m).Data(1:end-1)})
%                 else
%                     Li{2}=scatter(nanmean(TM_2(dyn,range,m)),nanmean(OG_2(dyn,range,m)),100,"filled",'MarkerFaceColor', 'r')  ;
%                     text(nanmean(TM_2(dyn,range,m))+.02,nanmean(OG_2(dyn,range,m)),{labels(m).Data(1:end-1)})
%                 end
 axis equal
        if dyn==1
%             title('Reactive')
%             xx=-.5:0.1:2.1;
%             plot(xx,xx,'r')

            %             xlim([-.5 2])
            %             ylim([-.5 2])
            yline(0)
            xline(0)
        else
            title('Contextual')
            yline(0)
            xline(0)
            %             xlim([-1 1])
            %             ylim([-1 1])
        end
%         xlabel('Treadmill')
%         ylabel('Overground')
xlabel('Reactive');ylabel('Contextual')
        set(gcf,'color','w')
    end

end
%%
%%Time Courses
colors=[0 0.4470 0.7410;0.8500 0.3250 0.0980];
adaptation=0
if adaptation==1
    range=1:480;
else
    range=481:680;
end

for m=1:28

    temp=[];
    x=[];
    figure
    for dyn=1:2
        if dyn==1
            x=-squeeze(TM(:,dyn,m,:));
            y=-squeeze(OG(:,dyn,m,:));
        else
            x=squeeze(TM(:,dyn,m,:));
            y=squeeze(OG(:,dyn,m,:));
        end
        x_mean=nanmean(x(range,:),2);
        y_mean=nanmean(y(range,:),2);

        P_x = prctile(x(range,:)',[2.5 97.5],1);
        P_y = prctile(y(range,:)',[2.5 97.5],1);
        %     x = 1:numel(y);
        %         x=index{i};
        %         std_dev = nanstd(d{i},1);

        subplot(2,1,1)
        hold on
        curve1 = P_y(1,:)';
        curve2 = P_y(2,:)';
        y2 = [1:size(x_mean,1), fliplr(1:size(x_mean,1))];
        inBetween = [curve1', fliplr(curve2')];
        fill(y2, inBetween,colors(dyn,:),'FaceAlpha',0.3,'EdgeColor','none')
        hold on;
        ylabel('W')
        xlabel('Strides')
        Li{dyn}=plot(1:size(x_mean,1),y_mean,'LineWidth', 2,'Color',colors(dyn,:));
        title(['Overground' ,{labels(m).Data(1:end-1)}])
        yline(0)

        subplot(2,1,2)
        hold on
        curve1 = P_x(1,:)';
        curve2 = P_x(2,:)';
        y2 = [1:size(x_mean,1), fliplr(1:size(x_mean,1))];
        inBetween = [curve1', fliplr(curve2')];
        fill(y2, inBetween,colors(dyn,:),'FaceAlpha',0.3,'EdgeColor','none')
        hold on;
        ylabel('W')
        xlabel('Strides')
        title(['Treadmill' ,{labels(m).Data(1:end-1)}])
        Li{dyn}=plot(1:size(x_mean,1),x_mean, 'LineWidth', 2,'Color',colors(dyn,:));
        yline(0)
        axis tight
    end

    legend([Li{:}],[{'Reactive';'Contextual'}])
    set(gcf,'color','w')
end
