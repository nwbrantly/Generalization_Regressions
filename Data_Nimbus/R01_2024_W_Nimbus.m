    
steps=41:45;

groupID={'CTR','NTR','VATR','CTS','NTS','VATS'};

figure
hold on
reactive=[];
contextual=[];

colors= [[0.8500 0.3250 0.0980];[0 0.4470 0.7410];[0.4660 0.6740 0.1880];[1 1 1];[1 1 1];[1 1 1]];
colorEdge= [[1 1 1];[1 1 1];[1 1 1];[0.8500 0.3250 0.0980];[0 0.4470 0.7410];[0.4660 0.6740 0.1880]];

for g=1:6
    
    load([groupID{g},'_post-adaptation_22-February-2024.mat'])
    reconstruction_contextual=squeeze(mean(contextual_trace(steps,:,:),'omitnan')); %contextual
    reconstruction_reactive=squeeze(mean(reactive_trace(steps,:,:),'omitnan')); % reactive
    
    for s=1:size(reconstruction_contextual,2)
        reconstruction_reactive(find(VIF_F(:,s)>5),s)=nan;
        reconstruction_reactive(find(R2{s}(:,1)<0.2),s)=nan;
        reconstruction_contextual(find(VIF_F(:,s)>5),s)=nan;
        reconstruction_contextual(find(R2{s}(:,1)<0.2),s)=nan;
        reconstruction_reactive(find(reconstruction_reactive(:,s)<0),s)=nan;
        reconstruction_contextual(find(reconstruction_contextual(:,s)<0),s)=nan;
        
    end
    
    reactive(:,g)= reshape(reconstruction_reactive,28*5,1);
    contextual(:,g)= reshape(reconstruction_contextual,28*5,1);
    n_muscle(g)=sum(~isnan(reactive(:,g)));
    n_muscle_C(g)=sum(~isnan(contextual(:,g)));
    subplot 211;hold on
    scatter(g,reactive(:,g),'o','filled','MarkerFaceColor',[.7 .7 .7])
    plot(g+.1,mean(reactive(:,g),'omitnan'),'sk','MarkerFaceColor',colors(g,:),...
        'MarkerSize',15,'MarkerEdgeColor',colorEdge(g,:))
    errorbar(g+.1,mean(reactive(:,g),'omitnan'),std(reactive(:,g),'omitnan'),'k')
    h= ttest(reactive(:,g));
    
    
    if h==1
        plot(g, 2, '*m')
    end
    
    subplot 212 ;hold on
    scatter(g,contextual(:,g),'o','filled','MarkerFaceColor',[.7 .7 .7])
    plot(g+.1,mean(contextual(:,g),'omitnan'),'sk','MarkerFaceColor',colors(g,:),...
        'MarkerSize',15,'MarkerEdgeColor',colorEdge(g,:))
    errorbar(g+.1,mean(contextual(:,g),'omitnan'),std(contextual(:,g),'omitnan'),'k')
    h= ttest(contextual(:,g));
    
    if h==1
        plot(g, 2, '*m')
    end
    
    
end

subplot 211;hold on
set(gca,'XTick',1:6,'XTickLabel',groupID,'FontSize',10)
yline(0)
ylabel({'Reactive';'(A.U)'})

subplot 212 ;hold on
set(gca,'XTick',1:6,'XTickLabel',groupID,'FontSize',10)
yline(0)
ylabel({'Contextual';'(A.U)'})
set(gcf,'color','w')

%%
%number of muscle per group
figure
hold on

for g=1:6
    
    bar(g,n_muscle(g),'FaceColor',colors(g,:),'EdgeColor',colorEdge(g,:))
    
end
title('Total of muscle per group')
ylabel('Number of muscles')
set(gca,'XTick',1:6,'XTickLabel',groupID,'FontSize',10)




%% Number of muscle with reactive and contextual > 0
figure
hold on

for g=1:6
    subplot 211; hold on
    bar(g,sum(reactive(:,g)>0),'FaceColor',colors(g,:),'EdgeColor',colorEdge(g,:))
    
    subplot 212; hold on
    bar(g,sum(contextual(:,g)>0),'FaceColor',colors(g,:),'EdgeColor',colorEdge(g,:))
    
end


subplot 211;hold on
title('W_{reactive} > 0')
set(gca,'XTick',1:6,'XTickLabel',groupID,'FontSize',10)
ylabel('Number of muscles')
yline(0)


subplot 212 ;hold on
title('W_{contextual} > 0')
set(gca,'XTick',1:6,'XTickLabel',groupID,'FontSize',10)
yline(0)
ylabel('Number of muscles')
set(gcf,'color','w')



%% ttest between conditions

[h_reactive,p_reactive] = ttest2(reactive(:,4),reactive(:,5))
[h_contextual,p_contextual] = ttest2(contextual(:,4),contextual(:,5))

%% Two-way ANOVA
CT=[contextual(:,1);contextual(:,4)];
NT=[contextual(:,2);contextual(:,5)];
VAT=[contextual(:,3);contextual(:,6)];

g1 = [ones(1,28*5) 2*ones(1,28*5) 3*ones(1,28*5) 1*ones(1,28*5) 2*ones(1,28*5) 3*ones(1,28*5) ];
g2=[ones(1,28*5*3) 2*ones(1,28*5*3)];

% reactive
c=reshape(reactive,28*5*6,1) ;
[~,~,stats] = anovan(c,{g1 g2},"Model","interaction", ...
    "Varnames",["Adaptation","PostAdaptation"]);

%contextual
c=reshape(contextual,28*5*6,1) ;
[~,~,stats] = anovan(c,{g1 g2},"Model","interaction", ...
    "Varnames",["Adaptation","PostAdaptation"]);
%% Stroke Data
%%

steps=41:45;

groupID={'C3TR_sub1','MWSTR_sub1','C3TS_sub1','MWSTS_sub1'};
figure
colors= [[1 1 1];[1 1 1];[0.8500 0.3250 0.0980];[0 0.4470 0.7410]];
colorEdge= [[0.8500 0.3250 0.0980];;[0 0.4470 0.7410];[1 1 1];[1 1 1];];


figure
hold on
reactive=[];
contextual=[];

% colors= [[1 1 1];[0 0.4470 0.7410];[0.8500 0.3250 0.0980]];
% colorEdge= [[0.8500 0.3250 0.0980];[1 1 1];[1 1 1];];

for g=1:length(groupID)
    
    load([groupID{g},'_post-adaptation_23-February-2024.mat'])
    reconstruction_contextual=squeeze(mean(contextual_trace(steps,:,:),'omitnan'))'; %contextual
    reconstruction_reactive=squeeze(mean(reactive_trace(steps,:,:),'omitnan'))'; % reactive
    
    for s=1%:size(reconstruction_contextual,2)
        reconstruction_reactive(find(VIF_F(:,s)>5),s)=nan;
        reconstruction_reactive(find(R2{s}(:,1)<0.2),s)=nan;
        reconstruction_contextual(find(VIF_F(:,s)>5),s)=nan;
        reconstruction_contextual(find(R2{s}(:,1)<0.2),s)=nan;
    end
    
    reactive(:,g)= reshape(reconstruction_reactive,28*1,1);
    contextual(:,g)= reshape(reconstruction_contextual,28*1,1);
    
    subplot 211;hold on
    scatter(g,reactive(:,g),'o','filled','MarkerFaceColor',[.7 .7 .7])
    plot(g+.1,mean(reactive(:,g),'omitnan'),'sk','MarkerFaceColor',colors(g,:),...
        'MarkerSize',15,'MarkerEdgeColor',colorEdge(g,:))
    errorbar(g+.1,mean(reactive(:,g),'omitnan'),std(reactive(:,g),'omitnan'),'k')
    h= ttest(reactive(:,g));
    
    
    if h==1
        plot(g, 2, '*m')
    end
    
    subplot 212 ;hold on
    scatter(g,contextual(:,g),'o','filled','MarkerFaceColor',[.7 .7 .7])
    plot(g+.1,mean(contextual(:,g),'omitnan'),'sk','MarkerFaceColor',colors(g,:),...
        'MarkerSize',15,'MarkerEdgeColor',colorEdge(g,:))
    errorbar(g+.1,mean(contextual(:,g),'omitnan'),std(contextual(:,g),'omitnan'),'k')
    h= ttest(contextual(:,g));
    
    if h==1
        plot(g, 2, '*m')
    end
    
    
end

subplot 211;hold on
set(gca,'XTick',1:6,'XTickLabel',groupID,'FontSize',10)
yline(0)
ylabel({'Reactive';'(A.U)'})

subplot 212 ;hold on
set(gca,'XTick',1:6,'XTickLabel',groupID,'FontSize',10)
yline(0)
ylabel({'Contextual';'(A.U)'})
set(gcf,'color','w')


g1 = [ones(1,28*1) 2*ones(1,28*1) 1*ones(1,28*1) 2*ones(1,28*1)];
g2=[ones(1,28*1*2) 2*ones(1,28*1*2)];

% reactive
c=reshape(reactive,28*1*4,1) ;
[~,~,stats] = anovan(c,{g1 g2},"Model","interaction", ...
    "Varnames",["Adaptation","PostAdaptation"]);

%contextual
c=reshape(contextual,28*1*4,1) ;
[~,~,stats] = anovan(c,{g1 g2},"Model","interaction", ...
    "Varnames",["Adaptation","PostAdaptation"]);

%%

steps=41:45;

groupID={'C3TR_sub9','C3TS_sub9'};

figure
hold on
reactive=[];
contextual=[];

colors= [[1 1 1];[0.8500 0.3250 0.0980]];
colorEdge= [[0.8500 0.3250 0.0980];[1 1 1]];

for g=1:length(groupID)
    
    load([groupID{g},'_post-adaptation_23-February-2024.mat'])
    reconstruction_contextual=squeeze(mean(contextual_trace(steps,:,:),'omitnan')); %contextual
    reconstruction_reactive=squeeze(mean(reactive_trace(steps,:,:),'omitnan')); % reactive
    
    for s=1:size(reconstruction_contextual,2)
        reconstruction_reactive(find(VIF_F(:,s)>5),s)=nan;
        reconstruction_reactive(find(R2{s}(:,1)<0.2),s)=nan;
        reconstruction_contextual(find(VIF_F(:,s)>5),s)=nan;
        reconstruction_contextual(find(R2{s}(:,1)<0.2),s)=nan;
    end
    
    reactive(:,g)= reshape(reconstruction_reactive,28*9,1);
    contextual(:,g)= reshape(reconstruction_contextual,28*9,1);
    
    subplot 211;hold on
    scatter(g,reactive(:,g),'o','filled','MarkerFaceColor',[.7 .7 .7])
    plot(g+.1,mean(reactive(:,g),'omitnan'),'sk','MarkerFaceColor',colors(g,:),...
        'MarkerSize',15,'MarkerEdgeColor',colorEdge(g,:))
    errorbar(g+.1,mean(reactive(:,g),'omitnan'),std(reactive(:,g),'omitnan'),'k')
    h= ttest(reactive(:,g));
    
    
    if h==1
        plot(g, 2, '*m')
    end
    
    subplot 212 ;hold on
    scatter(g,contextual(:,g),'o','filled','MarkerFaceColor',[.7 .7 .7])
    plot(g+.1,mean(contextual(:,g),'omitnan'),'sk','MarkerFaceColor',colors(g,:),...
        'MarkerSize',15,'MarkerEdgeColor',colorEdge(g,:))
    errorbar(g+.1,mean(contextual(:,g),'omitnan'),std(contextual(:,g),'omitnan'),'k')
    h= ttest(contextual(:,g));
    
    if h==1
        plot(g, 2, '*m')
    end
    
    
end

[h_reactive,p_reactive] = ttest2(reactive(:,1),reactive(:,2))
[h_contextual,p_contextual] = ttest2(contextual(:,1),contextual(:,2))

subplot 211;hold on
set(gca,'XTick',1:6,'XTickLabel',groupID,'FontSize',10)
yline(0)
ylabel({'Reactive';'(A.U)'})

subplot 212 ;hold on
set(gca,'XTick',1:6,'XTickLabel',groupID,'FontSize',10)
yline(0)
ylabel({'Contextual';'(A.U)'})
set(gcf,'color','w')
%%
muscleOrder={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT',...
    'VM', 'VL', 'RF', 'HIP','TFL', 'GLU'};
ytl=([strcat('f',muscleOrder) strcat('s',muscleOrder)]);  %List of muscle by fast and slow
ytl(end:-1:1) = ytl(:); % We flip the names to match the way labtools grab the data
yt=1:length(ytl); % getting the length of the vector
fs=14;
data=nan(length(ytl),2);
late=nan(length(ytl),2);

marker={'o','*','s'};
for t=1:3
    
    groupID={{'NTR','NTS'},{'CTR','CTS'},{'VATR','VATS'}};
    figure
    
    for g=1:2
        
        if g==1
            load([groupID{t}{1},'_post-adaptation_22-February-2024.mat'])
            x=1:2:56;
        else
            load([groupID{t}{2},'_post-adaptation_22-February-2024.mat'])
            x=1.5:2:57;
        end
        
        for s=1:5
            
            
            for m=1:28
                if R2{s}(m,1)>0.2 &&  VIF_F(m,s)<5
                    reactive(m,s)=mdl{m,1,s}.Coefficients.Estimate(1);
                    contextual(m,s)=mdl{m,1,s}.Coefficients.Estimate(2);
                else
                    reactive(m,s)=nan;
                    contextual(m,s)=nan;
                end
                
            end
            
            hold on
            
            if g==1
                subplot 211; hold on
                Li{1}=scatter(x, reactive(:,s), 'o','MarkerEdgeColor',[.7 .7 .7]);
                
                subplot 212; hold on
                Li{1}=scatter(x, contextual(:,s), 'o','MarkerEdgeColor',[.7 .7 .7]);
            else
                subplot 211; hold on
                Li{2}=scatter(x, reactive(:,s), 'o','filled','MarkerFaceColor',[.7 .7 .7]);
                subplot 212; hold on
                Li{2}=scatter(x, contextual(:,s), 'o','filled','MarkerFaceColor',[.7 .7 .7]);
            end
        end
        subplot 211; hold on
        %     scatter(x, nanmean(reactive,2), 'o');
        errorbar(x,nanmean(reactive,2),std(reactive,0,2,'omitnan'),"ok")
        temp(:,g)=nanmean(reactive,2);
        subplot 212; hold on
        errorbar(x,nanmean(contextual,2),std(contextual,0,2,'omitnan'),"ok")
        temp2(:,g)=nanmean(contextual,2);
    end
    subplot 211; hold on
    set(gca,'XTick',1:2:56,'XTickLabel',ytl,'FontSize',10)
    legend([Li{:}],{groupID{t}{1},groupID{t}{2}})
    yline(0)
    ylabel('W_{reactive}')
    subplot 212; hold on
    set(gca,'XTick',1:2:56,'XTickLabel',ytl,'FontSize',10)
    ylabel('W_{contextual}')
    yline(0)
    legend([Li{:}],{groupID{t}{1},groupID{t}{2}})
    
    set(gcf,'color','w')
    
    figure(10);
    hold on
    subplot 211;hold on
    Ll{t}=scatter(temp(:,1),temp(:,2),'Marker',marker{t});
    
    subplot 212; hold on
    Ll{t}=scatter(temp2(:,1),temp2(:,2),'Marker',marker{t});
    
    
end
x=-2:1:2;

figure(10);

hold on
subplot 211;hold on
plot(x,x,'k')
xlabel('W_{W Device - reactive}')
ylabel('W_{W/O Device - reactive}')
legend([Ll{:}],{'High','Low','Very Low'})

yline(0)
xline(0)


subplot 212;hold on
plot(x,x,'k')
xlabel('W_{W Device - contextual}')
ylabel('W_{W/O Device - contextual}')
yline(0)
xline(0)
set(gcf,'color','w')

%%


groupID={'C3TR_sub1','MWSTR_sub1','C3TS_sub1','MWSTS_sub1'};
figure
colors= [[1 1 1];[1 1 1];[0.8500 0.3250 0.0980];[0 0.4470 0.7410]];
colorEdge= [[0.8500 0.3250 0.0980];[0 0.4470 0.7410];[1 1 1];[1 1 1];];

for g=1:4
    load([groupID{g},'_post-adaptation_23-February-2024.mat'])
    
    for c=1:2
        n_muscle=0;
        for m=1:28
            data = mdl{m}.Coefficients.Estimate(c);
            SE =  mdl{m}.Coefficients.SE(c);
            
            upper = data + SE;
            lower = data - SE;
            
            
            if upper>0 && lower<0
                continue
            elseif  R2{1}(m,1)>0.2 &&  VIF_F(m,1)<5
                n_muscle = 1 + n_muscle;
            else
                continue
            end
            subplot(2,1,c)
            hold on
            bar(g,n_muscle,'FaceColor',colors(g,:),'EdgeColor',colorEdge(g,:))
            set(gca,'XTick',1:4,'XTickLabel',groupID,'FontSize',10)
            
        end
        
    end
end
subplot 211
ylabel('# of muscle reactive ')

subplot 212
ylabel('#  of muscle contextual')
set(gcf,'color','w')
%% Time courses and R2

binwith=5
analysis=0;isF=0;

% groupID={'VATR','VATS'};
% groupID={'NTR','NTS'};
%  groupID={'CTR','CTS'};
 groupID={'BAT'};
f=[1:14 1:14;1:14 1:14];

for g=1
%     load([groupID{g},'_post-adaptation_Indv_0_11-March-2024.mat'])
%  load([groupID{g},'_post-adaptation_Indv_0_13-March-2024.mat'])
  load([groupID{g},'_adaptation_Indv_0_14-March-2024.mat'])
    
    for muscle=1:28
        
        
        figure(f(g,muscle))
        
        % Pick the data that you want to plot
        Wdata(:,1)=[reactive_trace(41:150,muscle)]';
        Wdata(:,2)=[contextual_trace(41:150,muscle)]';
%         Wdata(:,3)=[baseline_trace(41:end,muscle)]';
        
        if muscle<15
            subplot(size(Wdata,2)+1,2,1)
        else
            subplot(size(Wdata,2)+1,2,2)
        end
        hold on
        scatter(1:length(movmean(Wdata(:,1),binwith)), movmean(Wdata(:,1),binwith),'filled','DisplayName',groupID{g}); %"#77AC30" )%
        
        legend
        
        title(labels(muscle).Data)
        
        yline(0,'HandleVisibility','off')
        ylabel({'Reactive';'(A.U)'})
        
        xlabel('strides')
        
        if size(Wdata,2)>=2
            if muscle<15
                subplot(size(Wdata,2)+1,2,3)
            else
                subplot(size(Wdata,2)+1,2,4)
            end
            hold on
            scatter(1:length(movmean(Wdata(:,2),binwith)), movmean(Wdata(:,2),binwith),'filled','DisplayName',groupID{g});
            legend
            yline(0,'HandleVisibility','off')
            ylabel({'Contextual';'(A.U)'})
            xlabel('strides')
        end
        
        if size(Wdata,2)>=3
            if muscle<15
                subplot(size(Wdata,2)+1,2,5)
            else
                subplot(size(Wdata,2)+1,2,6)
            end
            hold on
            scatter(1:length(movmean(Wdata(:,3),binwith)), movmean(Wdata(:,3),binwith),'filled','DisplayName',groupID{g});
            legend
            yline(0,'HandleVisibility','off')
            ylabel({'Baseline';'(A.U)'})
            xlabel('strides')
            
            
            if muscle<15
                subplot(size(Wdata,2)+1,2,7)
            else
                subplot(size(Wdata,2)+1,2,8)
            end

            
        else      
            
            if muscle<15
                subplot(size(Wdata,2)+1,2,5)
            else
                subplot(size(Wdata,2)+1,2,6)
            end
     
   
        end
        % Plot the time course of the R2
            r2 = my_Rsquared_coeff(model{muscle}.EMGobserved', model{muscle}.EMG_estimated,1);
            hold on
            plot(movmean(r2(41:150),binwith),'DisplayName',groupID{g})
            ylim([-.5 1])
            yline(0,'HandleVisibility','off')
            legend
            xlabel('strides')
            ylabel('R^2')
            set(gcf,'color','w')
        
        if muscle<=14
            isF=0;
        else
            isF=1;
        end
        
        
        
    end
    
end