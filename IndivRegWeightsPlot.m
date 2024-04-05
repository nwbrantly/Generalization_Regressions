close all; clc
Opacity = 0.3;
colorOrder = PlotHelper.colorOrder;

%% Time courses
labels = ytl;
% weights(:,1,:,:) = reactive_trace;
% weights(:,2,:,:) = contextual_trace;

weights(:,1,:,:) = dynamics_TR(:,1,:,:);
weights(:,2,:,:) = dynamics_TR(:,2,:,:);

E=std(weights,0,4,'omitnan');%./sqrt(size(weights,4)); %Standar error, std, 2nd arg specifies normalization by N-1
avgOverSub = nanmean(weights,4);
strides = 41:size(weights,1); %argument needs to be in row vector format: 1xstrides

for mIdx = 1:28
    if mod(mIdx,8) == 1
        f = figure('units','normalized','outerposition',[0 0 1 1])
        subPltIdx = 1;
    end
    subplot(4,2,subPltIdx);
    color = colorOrder(1,:);
    [Pa, Li]= nanJackKnife(strides,avgOverSub(strides,1,mIdx)',E(strides,1,mIdx)',color,color+0.5.*abs(color-1),Opacity);
    pp=patch([40 480 480 40],[-1 -1 1 1],.7*ones(1,3),'FaceAlpha',.5,'EdgeColor','none'); %PATR - PATS
    uistack(pp,'bottom')
    yline(0);
    set(gcf,'color','w')
    %     xline([40,930,1370]);
    title(['Reactive ' labels{mIdx}]);
    xlabel('Stride')
    
    
    subplot(4,2,subPltIdx);
    color = colorOrder(2,:);
    [Pa, Li]= nanJackKnife(strides,avgOverSub(strides,2,mIdx)',E(strides,2,mIdx)',color,color+0.5.*abs(color-1),Opacity);
    title(['Contextual ' labels{mIdx}]);
    xlabel('Stride')
    pp=patch([40 480 480 40],[-1 -1 1 1],.7*ones(1,3),'FaceAlpha',.5,'EdgeColor','none'); %PATR - PATS
    uistack(pp,'bottom')
    yline(0);
    set(gcf,'color','w')
    %     xline([40,930,1370]);
    %     sgtitle(labels{mIdx})
    subPltIdx = subPltIdx +1;
end

%% Time course two groups
muscleOrder={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT',...
    'VM', 'VL', 'RF', 'HIP','TFL', 'GLU'};
ytl=([strcat('f',muscleOrder) strcat('s',muscleOrder)]);  %List of muscle
ytl(end:-1:1) = ytl(:);
yt=1:length(ytl);

labels = ytl;

%%
%CTS group
load('CTS_Indv_1_regressorConstant_0_removeBias_1_PostAdaptationV2.mat')
%put adapt and OGPost together in order: OGBase (40), TMBase(40), Adapt(450), OGPost(300)
%(with a minus sign in front)
%adaptation take only the first 340 strides (40 baseline and 300
%adaptation) for the 2 groups to match
% ada900.yhat(:,1:340,:)
% group1.reactive_trace = cat(1, indivSubOGPost.reactive(1:40,:,:), indivSubAda.reactive(:,:,:), -indivSubOGPost.reactive(41:end,:,:));
% group1.contextual_trace = cat(1, indivSubOGPost.context(1:40,:,:), indivSubAda.context(:,:,:), indivSubOGPost.context(41:end,:,:));
group1.reactive_trace = [reactive_trace(1:40,:,:); -reactive_trace(40:end,:,:)];
group1.contextual_trace = [contextual_trace(1:40,:,:); contextual_trace(40:end,:,:)];


% NTS group
load('NTS_Indv_1_regressorConstant_0_removeBias_1_PostAdaptationV2.mat')
%put adapt and OGPost together in order: OGBase(40), TMBase(40),
%Adapt(900), OGPost (300)
%(with a minus sign in front)
% group2.reactive_trace = cat(1, indivSubOGPost.reactive(1:40,:,:), indivSubAda.reactive(:,:,:), -indivSubOGPost.reactive(41:end,:,:));
% group2.contextual_trace = cat(1, indivSubOGPost.context(1:40,:,:), indivSubAda.context(:,:,:), indivSubOGPost.context(41:end,:,:));
group2.reactive_trace = [reactive_trace(1:40,:,:); -reactive_trace(40:end,:,:)];
group2.contextual_trace = [contextual_trace(1:40,:,:); contextual_trace(40:end,:,:)];

groupNames = {'CTS','NTS'};

% plot time course
%set up plot configurations.
Opacity = 0.3;
colorOrder = PlotHelper.colorOrder;
moveStrides = 5;
% conditionSepStrides = [40,80,380]; %end of OGBase, end of TMBase, end of adaptation
% conditionSepString = {'TMBase','Ada','OGPost'};
conditionSepStrides = [0 41]; %new version plot with gap, start of ada, SS, postadap
conditionSepString = {'Base','OGPost'};
statesName = {'Reactive','Context'};

% weights = dataBySubjects.weights_tranp;
for groupIdx = 1:length(groupNames)
    weights = [];
    weights(:,1,:,:) = eval(['group' num2str(groupIdx) '.reactive_trace']);
    weights(:,2,:,:) =eval(['group' num2str(groupIdx) '.contextual_trace']);%  contextual_trace;
    E=movstd(std(weights,0,4,'omitnan')./sqrt(size(weights,4)),moveStrides);
    avgOverSub = movmean(nanmean(weights,4),moveStrides);
    eval(['group'  num2str(groupIdx) '.E=E;']) %Standar error, std, 2nd arg specifies normalization by N-1
    eval(['group' num2str(groupIdx) '.avgOverSub = avgOverSub;'])
end
% strides = 41:size(weights,1); %argument needs to be in row vector format: 1xstrides
strides = 1:size(weights,1); %argument needs to be in row vector format: 1xstrides

%AUF protocl: strides=[-40 890 440 140]; %Number per strides per condition % cond={'TMBase','Adaptation','OGPost', 'TMPost'}; %Conditions for this group
for mIdx = 1:28
    if mod(mIdx,6) == 1
        f = figure('units','normalized','outerposition',[0 0 1 1])
        subPltIdx = 1;
    end
    for stateIdx = 1:2 %reactive, then context
        subplot(3,4,subPltIdx+stateIdx-1); hold on;
        for groupIdx = 1:length(groupNames)
            color = colorOrder(groupIdx,:);
            avgOverSub = eval(['group' num2str(groupIdx) '.avgOverSub']);
            E = eval(['group' num2str(groupIdx) '.E']);
            %TODO: plot in 2 separate batches: ada, then steady states ada, and then OGPost with a gap
            %in between the strides;
            %TODO: handle different strides of the 2 groups.
            [Pa, Li{groupIdx}]= nanJackKnife(1:40,avgOverSub(1:40,stateIdx,mIdx)',E(1:40,stateIdx,mIdx)',color,color+0.5.*abs(color-1),Opacity);
            %last 40 of ada, xindex shift by 50 to have a gap
            %             [Pa, Li{groupIdx}]= nanJackKnife(41:341,avgOverSub(end-339:end-300,stateIdx,mIdx)',E(end-339:end-300,stateIdx,mIdx)',color,color+0.5.*abs(color-1),Opacity);
            %300 strides of post adaptation, xindex shift by 50 again for a gap
            [Pa, Li{groupIdx}]= nanJackKnife(41:340,avgOverSub(end-299:end,stateIdx,mIdx)',E(end-299:end,stateIdx,mIdx)',color,color+0.5.*abs(color-1),Opacity);
        end
        %     color = colorOrder(3,:);
        %     [Pa, Li2]= nanJackKnife(strides,group2.avgOverSub(strides,1,mIdx)',group2.E(strides,1,mIdx)',color,color+0.5.*abs(color-1),Opacity);
        xlim([-5 350])
        yline(0);
        xline(conditionSepStrides,'-',conditionSepString,'LabelHorizontalAlignment','center'); %,'LabelVerticalAlignment','bottom'
        title([statesName{stateIdx} ' ' labels{mIdx}]);
        xlabel('Stride')
        if subPltIdx == 1
            legend([Li{:}],groupNames);
        end
    end
    
    %     subplot(3,4,subPltIdx+1); hold on;
    %     for groupIdx = 1:length(groupNames)
    %         color = colorOrder(groupIdx,:);
    %         avgOverSub = eval(['group' num2str(groupIdx) '.avgOverSub']);
    %         E = eval(['group' num2str(groupIdx) '.E']);
    %         [Pa, Li{groupIdx}]= nanJackKnife(strides,avgOverSub(strides,2,mIdx)',E(strides,2,mIdx)',color,color+0.5.*abs(color-1),Opacity);
    %     end
    % %     color = colorOrder(1,:);
    % %     [Pa, Li1]= nanJackKnife(strides,group1.avgOverSub(strides,2,mIdx)',group1.E(strides,2,mIdx)',color,color+0.5.*abs(color-1),Opacity);
    % %     color = colorOrder(3,:);
    % %     [Pa, Li2]= nanJackKnife(strides,group2.avgOverSub(strides,2,mIdx)',group2.E(strides,2,mIdx)',color,color+0.5.*abs(color-1),Opacity);
    %
    %     title(['Contextual ' labels{mIdx}]);
    %     xlabel('Stride')
    %     yline(0);
    %     xline(conditionSepStrides,'-',{'TMBase','Ada','OGPost'},'LabelHorizontalAlignment','center'); %,'LabelVerticalAlignment','bottom'
    %     if subPltIdx == 1 %add legend in the end aftre all graphical component are added, skip vertical lines legend.
    %         legend([Li{:}],groupNames);
    %     end
    
    subPltIdx = subPltIdx +2;
    set(gcf,'color','w')
end

%% plot r2 barplots

figure
hold on

names = fieldnames(r2);

for type=1:length(fieldnames(r2))
    
    hold on
    bar(type, nanmean(r2.(names{type})))
    %     errorbar(type, nanmean(r2.(names{type})), std(r2.(names{type}))/sqrt(length(r2.(names{type}))),'k')
    errorbar(type, nanmean(r2.(names{type})), std(r2.(names{type})),'k')
    
    % for s=1:length( r2.(names{type}))
    
    scatter(type+.1, r2.(names{type}),25,"filled",'MarkerFaceColor', 'k')
    
    % end
end

set(gca,'XTick', 1:length(fieldnames(r2)) ,'XTickLabel',names);
ylabel('R^2')
title('Post-Adapt')
set(gcf,'color','w')

%% plot VAF

figure
hold on
li=[];
names = fieldnames(r2uncenter);


load('CTS_Indv_1_regressorConstant_0_removeBias_1_PostAdaptationV2.mat')
r2uncenter2=r2uncenter;

load 'NTS_Indv_1_regressorConstant_0_removeBias_1_PostAdaptationV2.mat'

for type=1:length(fieldnames(r2uncenter))
    type2=1:2:6; %1:3;%
    
    hold on
    
    %%% My data
    li{1}=bar(type2(type), nanmean(r2uncenter.(names{type})),'FaceColor',"#D95319");
    %     errorbar(type, nanmean(r2uncenter.(names{type})), std(r2uncenter.(names{type}))/sqrt(length(r2uncenter.(names{type}))),'k')
    errorbar(type2(type), nanmean(r2uncenter.(names{type})), std(r2uncenter.(names{type})),'k')
    
    % for s=1:length(r2uncenter.(names{type}))
    
    scatter(type2(type)+.1, r2uncenter.(names{type}),25,"filled",'MarkerFaceColor', 'k')
    
    % end
    
    %% Shuqi Data
    shuqi=2:2:6;
    %     li{2}=bar(shuqi(type), nanmean(R2.shifted1.relativeToMean0(type,:)),'FaceColor',"#D95319");
    %     errorbar(shuqi(type), nanmean(R2.shifted1.relativeToMean0(type,:)), nanstd(R2.shifted1.relativeToMean0(type,:)),'k')
    %     scatter(shuqi(type)+.1, R2.shifted1.relativeToMean0(type,:),25,"filled",'MarkerFaceColor', 'k')
    li{2}=bar(shuqi(type), nanmean(r2uncenter2.(names{type})),'FaceColor',"#0072BD");
    errorbar(shuqi(type), nanmean(r2uncenter2.(names{type})), std(r2uncenter2.(names{type})),'k')
    scatter(shuqi(type)+.1, r2uncenter2.(names{type}),25,"filled",'MarkerFaceColor', 'k')
    
end


legend([li{:}],{'NTS','CTS'})
set(gca,'XTick', 1.5:2:length(fieldnames(r2))*2 ,'XTickLabel',names);
ylabel('VAF')
title('Post-Adapt')
set(gcf,'color','w')



%% plot W's contextual
color = PlotHelper.colorOrder;

for g=1:2
    if g==1
        load('BATS_Indv_0_Adapt_0_regressorConstant_0_removeBias_1.mat')
    elseif g==2
        load('/Users/dulcemariscal/Downloads/group_model_noIntercept_OGPost.mat')
    end
    
    for i=1:28
        
        subplot 211
        hold on
        li{g}=scatter(i,mdl{i}.Coefficients.Estimate(1),50,"filled",'MarkerFaceColor', color(g,:));
        errorbar(i,mdl{i}.Coefficients.Estimate(1),mdl{i}.Coefficients.SE(1),'k')
        
        subplot 212
        hold on
        li2{g}=scatter(i,mdl{i}.Coefficients.Estimate(2),50,"filled",'MarkerFaceColor', color(g,:));
        errorbar(i,mdl{i}.Coefficients.Estimate(2),mdl{i}.Coefficients.SE(2),'k')
        
    end
    
end

subID={'YA','OA'};
subplot 211
ylabel({'Reactive';'AU'})
yline(0)
set(gca,'XTick',yt,'XTickLabel',ytl,'FontSize',10)
legend([li{:}],subID, 'location','Best')
set(gcf,'color','w')

subplot 212
ylabel({'Contextual';'AU'})
yline(0)
set(gca,'XTick',yt,'XTickLabel',ytl,'FontSize',10)
legend([li2{:}],subID, 'location','Best')
set(gcf,'color','w')

%% W's bootstrapping

% load BAT_24_iteration_2000_Individual_muscles_C_constant_PostAdaptation_per_group_24-Jul-2023.mat
load('BAT_24_iteration_2000_Individual_muscles_C_constant_PostAdaptation_per_group_25-July-2023 08:24.mat')
dynamics_TS=dynamics_TS;
post1Index = 481:485;

load AUFV02_22_iteration_2000_Individual_muscles_C_constant_Adaptation_per_group.mat
dynamics_TS=dynamics_TR;
post1Index = 931:935;

muscleOrder={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT',...
    'VM', 'VL', 'RF', 'HIP','TFL', 'GLU'};
ytl=([strcat('f',muscleOrder) strcat('s',muscleOrder)]);  %List of muscle
ytl(end:-1:1) = ytl(:);
yt=1:length(ytl);


titles = {'Reactive','Contextual'};

for cIdx=1:2 %reactive, then context,
    
    f = figure('units','normalized','outerposition',[0 0 1 0.5]);%('Position', get(0, 'Screensize'));
    
    for mIdx = 1:size(dynamics_TS,3) %for each muscle  muscle orders start from sGLU
        postData = nanmean(squeeze(dynamics_TS(post1Index,cIdx,mIdx,:)));
        PlotHelper.plotCI(mIdx, postData,'b',ytl{mIdx},false, true) %
    end
    
    title([titles{cIdx} ' (Mean, 95% from bootstrap)'])
    xlim([0 size(dynamics_TR,3)+1])
    xticks(1:size(dynamics_TR,3))
    xticklabels(ytl)
    xtickangle(45)
    PlotHelper.tightMargin(gca)
    
    set(gcf,'color','w')
end


%% W's bootstrapping with older adults


f(1)=figure('units','normalized','outerposition',[0 0 1 0.5]);
f(2)=figure('units','normalized','outerposition',[0 0 1 0.5]);

muscleOrder={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT',...
    'VM', 'VL', 'RF', 'HIP','TFL', 'GLU'};
ytl=([strcat('f',muscleOrder) strcat('s',muscleOrder)]);  %List of muscle
ytl(end:-1:1) = ytl(:);
yt=1:length(ytl);

titles = {'Reactive','Contextual'};

for cIdx=1:2 %reactive, then context,
    
    for g=1:2 % young and then old adults
        
        if g==1
%             load BAT_24_iteration_2000_Individual_muscles_C_constant_PostAdaptation_per_group_24-Jul-2023.mat
load 'BAT_24_iteration_2000_Individual_muscles_C_constant_PostAdapt_per_group_19-December-2023.mat'
            %                 load('BAT_24_iteration_2000_Individual_muscles_C_constant_PostAdaptation_per_group_25-July-2023 08:24.mat')
            dynamics=dynamics_TS;
            post1Index = 41:45;
        else
%             load AUFV02_22_iteration_2000_Individual_muscles_C_constant_PostAdaptation_per_group.mat
%             dynamics_TS=dynamics_TR;
%             post1Index = 931:935;
                             dynamics=dynamics_TS;
                             post1Index = 220:240;
            
        end
        
        figure(cIdx) %= figure('units','normalized','outerposition',[0 0 1 0.5]);%('Position', get(0, 'Screensize'));
        hold on
        
        for mIdx = 1:size(dynamics,3) %for each muscle  muscle orders start from sGLU
            
            postData = nanmean(squeeze(dynamics(post1Index,cIdx,mIdx,:)));
            yMean = mean(postData);
            yCI95= prctile(postData,[2.5 97.5]);
            
            if g==1
                x=1:2:28*2;
                lineColor = 'b';
            else
                x=2:2:28*2;
                lineColor = 'k';
            end
            
            plot([x(mIdx),x(mIdx)], yCI95,'Color',lineColor,'LineWidth',3.5,'Marker','None','HandleVisibility','off')
            
            if (yCI95(1) > 0 || yCI95(2) < 0)
                %                 plot(x(mIdx),yCI95(2)+0.2,'*r');
                li{3}=plot(x(mIdx),yMean,'.','Color','r','LineWidth',3.5,'MarkerSize',38,'HandleVisibility','off');
            else
                li{g}=plot(x(mIdx),yMean,'.','Color',lineColor,'LineWidth',3.5,'MarkerSize',38,'HandleVisibility','off');
            end
            %              plot(x(mIdx),yMean,'.','Color',lineColor,'LineWidth',3.5,'MarkerSize',38,'HandleVisibility','off');
            
            
        end
        
        title([titles{cIdx} ' (Mean, 95% from bootstrap)'])
        xlim([0 size(dynamics,3)*2+1])
        xticks(1.5:2:size(dynamics,3)*2)
        xticklabels(ytl)
        xtickangle(45)
        %     PlotHelper.tightMargin(gca)
        yline(0)
        set(gcf,'color','w')
        
    end
    
    for x=2.5:2:28*2
        xline(x)
    end
    %       legend([li{:}],{'YA','OA','\neq0'},'Location','Best')
    legend([li{:}],{'Early','Late','\neq0'},'Location','Best')
    set(findall(gcf,'-property','FontSize'),'FontSize',25)
end

%      legend
