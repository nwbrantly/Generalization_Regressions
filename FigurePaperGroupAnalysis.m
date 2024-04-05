% Figure paper based on group analysis.
clear all
%% Organizaing muscle list. This is necesary for the labeling later if we remove muscle
muscleOrder={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT',...
    'VM', 'VL', 'RF', 'HIP','TFL', 'GLU'};
ytl=([strcat('f',muscleOrder) strcat('s',muscleOrder)]);  %List of muscle by fast and slow
ytl(end:-1:1) = ytl(:); % We flip the names to match the way labtools grab the data
yt=1:length(ytl); % getting the length of the vector
fs=14;

% Paired t-test and Wilcoxon Signed Rank Test for adaptation data (Early vs. Late)

load('BAT_adaptation_31-January-2024.mat') %This matrix needs to be uploaded to the server.
% load('BATS_post-adaptation_31-January-2024.mat')
yt=1:length(ytl);
early=nan(length(ytl),2);
late=nan(length(ytl),2);
figure 


 xtl={'Early_{adaptation}','Late_{adaptation}'};
for m= 1:length(yt)
    muscle= m; %Needed to skip muscles with VIF>5
    if le(0,R2{1}(muscle,1)) &&  le(0,R2{1}(muscle,2)) &&  VIF_F(muscle,1)<5 % Are early and late R2>0 and VIF<5?
        early(m,:)=mdl{muscle,1}.Coefficients.Estimate';
        
        late(m,:)=mdl{muscle,2}.Coefficients.Estimate';
        
        subplot 325
        hold on 
%         plot(m,R2{1}(m,1),'o', "MarkerFaceColor",[0.6350 0.0780 0.1840],"MarkerEdgeColor",[0.6350 0.0780 0.1840])
%         ylabel('R^2')
        yyaxis left
        scatter(m,R2{1}(m,1),'o','filled')
        ylabel('R^2')
        yyaxis right
        scatter(m,VIF_F(m,1),'o','filled')
%         plotyy(m,R2{1}(m,2),m,VIF_F(m,1))
        ylabel('VIF')
        yline(0)
        title(xtl{1})
        set(gca,'XTick',yt,'XTickLabel',ytl,'FontSize',10)
        
        subplot 326
        hold on
        yyaxis left
        scatter(m,R2{1}(m,2),'o','filled')
        ylabel('R^2')
        yyaxis right
        scatter(m,VIF_F(m,1),'o','filled')
%         plotyy(m,R2{1}(m,2),m,VIF_F(m,1))
        ylabel('VIF')
        yline(0)
        title(xtl{2})
        set(gca,'XTick',yt,'XTickLabel',ytl,'FontSize',10)
    end
end


[~,p(1),ci{1},stats{1}]  = ttest(early(:,1),late(:,1),'Tail','right');
[pnp_reactive,~,statsnp_reactive]  = signrank(early(:,1),late(:,1),'Tail','right');



[~,p(2),ci{2},stats{2}]  = ttest(early(:,2),late(:,2),'Tail','left');
[pnp_context,~,statsnp_context]  = signrank(early(:,2),late(:,2),'Tail','left');


% Adaptation Early vs Late
% close all
% load('BAT_adaptation_31-January-2024.mat')
% Reactive
cond={'Reactive','Contextual'};
for c=1:2
    
    subplot(3,1,c)
    early_mean= mean(early(:,c),'omitnan');
    late_mean= mean(late(:,c),'omitnan');
    early_std= std(early(:,c),'omitnan');
    late_std= std(late(:,c),'omitnan');
   
    %     figure
    hold on
    for m=1:length(early)
        
        plot([1 2],[early(m,c) late(m,c)],'-o','Color',[.7 .7 .7])
        
    end
    plot([1 2],[early_mean late_mean],'-ok')
%     errorbar([early_mean late_mean], [early_std late_std ]/sqrt(length(early(~isnan(early(:,1))))),'k')
%     
    
    ylabel([cond{c}])
    xlim([.75 2.25])
    xticks([1 2])
    xticklabels(xtl)
    
    %
    yline(0)
    % legend(muscleOrder(end:-1:1),'Location','Best')
    % set(findall(gcf,'-property','FontSize'),'FontSize',25)
    set(gcf,'color','w')
    %      PlotHelper.tightMargin(gca)
%     if p(c)<0.05
%         
%         plot(1.5, max([early(:,c) late(:,c)],[],'all')+.1,'*r')
%         plot([1 2], [max([early(:,c) late(:,c)],[],'all')+0.05  max([early(:,c) late(:,c)],[],'all')+.05] ,'k','LineWidth',2)
%     end
    
end






%% Post-adaptation
clear all; close all
muscleOrder={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT',...
    'VM', 'VL', 'RF', 'HIP','TFL', 'GLU'};
ytl=([strcat('f',muscleOrder) strcat('s',muscleOrder)]);  %List of muscle by fast and slow
ytl(end:-1:1) = ytl(:); % We flip the names to match the way labtools grab the data
yt=1:length(ytl); % getting the length of the vector
fs=14;
data.groups{1}=nan(28,2);
data.groups{2}=nan(28,2);
early=nan(length(ytl),1);
groups={'TM','OG'};
figure
hold on
l=0;
plotOrder=[1,5,2,6];
xtl={'TM','OG'};

for g=1:2
    if g==1
        load('BATR_post-adaptation_31-January-2024.mat')
        VIF(:,1)=VIF_F;
    else
        load('BATS_post-adaptation_31-January-2024.mat')
        VIF(:,2)=VIF_F;
        
    end
    
    
    % reactive
    cond={'Reactive','Contextual'};
    for c=1:2
        l=l+1;
        subplot(5,2,plotOrder(l))
        hold on
        for m= 1:length(yt)
            muscle= m; %Needed to skip muscles with VIF>5
            
            if le(0,R2{1}(muscle,1))  &&  VIF_F(muscle,1)<5 % Are early and late R2>0 and VIF<5?
                data.groups{g}(m,c)=mdl{muscle,1}.Coefficients.Estimate(c)';
                errorbar(muscle, mdl{muscle,1}.Coefficients.Estimate(c),mdl{muscle,1}.Coefficients.SE(c),'k')
                if c==1
                    plot(muscle, data.groups{g}(m,c),'o', "MarkerFaceColor",[0.9290 0.6940 0.1250],"MarkerEdgeColor",'k');
                else
                    plot(muscle, data.groups{g}(m,c),'o', "MarkerFaceColor",[0 0.4470 0.7410],"MarkerEdgeColor",'k');
                end
                
            elseif  VIF_F(muscle,1)<5  % Are early and late R2>0 and VIF<5?
                errorbar(muscle, mdl{muscle,1}.Coefficients.Estimate(c),mdl{muscle,1}.Coefficients.SE(c),'k')
                if c==1
                    plot(muscle,mdl{muscle,1}.Coefficients.Estimate(c),'o', "MarkerFaceColor",[.7 .7 .7],"MarkerEdgeColor",'k');
                else
                    plot(muscle,mdl{muscle,1}.Coefficients.Estimate(c),'o', "MarkerFaceColor",[.7 .7 .7],"MarkerEdgeColor",'k');
                end
                
                
            end
        end
        ylabel([cond{c}])
        yline(0)
        title(groups{g})
        
        if c==1
            ylim([0 1.6])
        else
            ylim([-.5 1])
        end
        set(gca,'XTick',yt,'XTickLabel',ytl,'FontSize',10)
    end
    
    subplot(5,2,8+g)
    plot(1:length(R2{1}),R2{1}(:,1),'o', "MarkerFaceColor",[0.6350 0.0780 0.1840],"MarkerEdgeColor",[0.6350 0.0780 0.1840])
    ylabel('R^2')
    yline(0)
    title(groups{g})
    set(gca,'XTick',yt,'XTickLabel',ytl,'FontSize',10)
    
end

temp=[2,4];
p=[];

for c=1:2
    TM=[];OG=[];
    subplot(5,1,temp(c))
    hold on
    for  m= 1:length(yt)
        if ~isnan(data.groups{1}(m,c)) && ~isnan(data.groups{2}(m,c)) && VIF(m,1)<5 && VIF(m,2)<5  % Are there NAN in the TM or OG data and VIF<5?
            plot([1 2], [data.groups{1}(m,c) data.groups{2}(m,c)],'-o','Color',[.7 .7 .7])
            TM=[TM; data.groups{1}(m,c)];
            OG=[OG; data.groups{2}(m,c)];
            
            xlim([.75 2.25])
            ylabel([cond{c}])
            yline(0)
            
            if c==1
                ylim([0 1.5])
            else
                ylim([-.2 1])
            end
            
        end
    end
    plot([1 2],[mean(TM,'omitnan') mean(OG,'omitnan')],'-ok')
    errorbar([mean(TM,'omitnan') mean(OG,'omitnan')], [std(TM,'omitnan') std(OG,'omitnan')]/sqrt(length(early(~isnan(OG)))),'k')
    
    
    xlim([.9 2.1])
    xticks([1 2])
    xticklabels(xtl)
    [~,p]  = ttest2(TM,OG)
    %      [p_np]  = signrank(TM,OG,'Tail','right')
    
    if p<0.05
        
        plot(1.5, max([data.groups{1}(:,c),data.groups{2}(:,c)],[],'all')+.1,'*r')
        plot([1 2], [max([TM OG],[],'all')+0.05 max([TM OG] ,[],'all')+.05],'k','LineWidth',2)
    end
    
end


set(gcf,'color','w')

% [~,p(1),ci{1},stats{1}]  = ttest(data.groups{1}(:,1),data.groups{2}(:,1),'Tail','right')
% [pnp_reactive,~,statsnp_reactive]  = signrank(data.groups{1}(:,1),data.groups{2}(:,1),'Tail','right')
%
% [~,p(2),ci{2},stats{2}]  = ttest(data.groups{1}(:,2),data.groups{2}(:,2),'Tail','left')
% [pnp_context,~,statsnp_context]  = signrank(data.groups{1}(:,2),data.groups{2}(:,2),'Tail','left')

%% Scatter plots


clear all; close all
muscleOrder={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT',...
    'VM', 'VL', 'RF', 'HIP','TFL', 'GLU'};
ytl=([strcat('f',muscleOrder) strcat('s',muscleOrder)]);  %List of muscle by fast and slow
ytl(end:-1:1) = ytl(:); % We flip the names to match the way labtools grab the data
yt=1:length(ytl); % getting the length of the vector
fs=14;
data.groups{1}=nan(28,2);
data.groups{2}=nan(28,2);
early=nan(length(ytl),1);
groups={'TM','OG'};

l=0;
plotOrder=[1,5,2,6];
xtl={'TM','OG'};



load('BATR_post-adaptation_31-January-2024.mat')
VIF=VIF_F;
data_TM= mdl(:,1);
R2_TM=R2{1}(:,1);

load('BATS_post-adaptation_31-January-2024.mat')
data_OG= mdl(:,1);
R2_OG=R2{1}(:,1);

hip=[0.4660, 0.6740, 0.1880];
thigh=   [0.4940, 0.1840, 0.5560];
shank= 	[0.8500, 0.3250, 0.0980];
colors=[hip; hip;hip; thigh; thigh; thigh; thigh; thigh; thigh; shank; shank; shank; shank; shank ] ;
colors=[colors;colors];
plotOrder=[1,2;5,6];

cond={'Reactive','Contextual'};
l=0;

for c=1:2
    %      figure
    l=l+1;
    %       subplot(4,2,plotOrder(l))
    for m= 1:length(yt)
        
        if m<14
            subplot(4,2, plotOrder(c,1));
            hold on
            
            
            
            
        else
            subplot(4,2, plotOrder(c,2));
            hold on
            
            
            
        end
        
        if VIF(m,1)<5
            
            if le(0, R2_OG(m,1)) && le(0, R2_TM(m,1)) %If the TM and OG R2>0
                
                %                  scatter( data_TM{m}.Coefficients.Estimate(c)', data_OG{m}.Coefficients.Estimate(c)')
                errorbar(data_TM{m}.Coefficients.Estimate(c)', data_OG{m}.Coefficients.Estimate(c)',...
                    data_TM{m}.Coefficients.SE(c),data_TM{m}.Coefficients.SE(c),data_OG{m}.Coefficients.SE(c),...
                    data_OG{m}.Coefficients.SE(c),'o-','Color', colors(m,:))
                
            else
                %                   scatter( data_TM{m}.Coefficients.Estimate(c)', data_OG{m}.Coefficients.Estimate(c)','filled','MarkerFaceColor',[0.7 .7 .7])
                errorbar(data_TM{m}.Coefficients.Estimate(c)', data_OG{m}.Coefficients.Estimate(c)',...'
                    data_TM{m}.Coefficients.SE(c),data_TM{m}.Coefficients.SE(c),data_OG{m}.Coefficients.SE(c),...
                    data_OG{m}.Coefficients.SE(c),...
                    'o-', 'Color',[.7 .7 .7])
            end
            
        end
        
        
    end
    
    
end
for c=1:2
    for f=1:2
        if f==1
            subplot(4,2, plotOrder(c,f));
            
            
            if c==1
                y=-.5:.1:1.5;
                axis([-.3 1.5 -.3 1.5])
            else
                y=-.5:.1:1;
                axis([-.3 1 -.3 1])
            end
            
            plot(y,y,'k')
            
            title(['s' ,cond{c}])
            xlabel('W_{TM}');ylabel('W_{OG}')
            set(gcf,'color','w')
            yline(0,'-.k');xline(0,'-.k')
        else
            subplot(4,2, plotOrder(c,f));
            plot(y,y,'k')
            title(['f' ,cond{c}])
            xlabel('W_{TM}');ylabel('W_{OG}')
            set(gcf,'color','w')
            
            if c==1
                y=-.5:.1:1.5;
                axis([-.3 1.5 -.3 1.5])
            else
                y=-.5:.1:1;
                axis([-.3 1 -.3 1])
            end
            yline(0,'-.k');xline(0,'-.k')
        end
        
        
        
        
    end
end

temp=[2,4];

p=[];
for c=1:2
    TM=[];OG=[];
    subplot(4,1,temp(c))
    hold on
    for  m= 1:length(yt)
        if  VIF(m,1)<5 &&  R2_OG(m,1)>0 && R2_TM(m,1)>0 % Are there NAN in the TM or OG data and VIF<5?
            plot([1 2], [data_TM{m}.Coefficients.Estimate(c) data_OG{m}.Coefficients.Estimate(c)],'-o','Color',[.7 .7 .7])
            TM=[TM; data_TM{m}.Coefficients.Estimate(c)];
            OG=[OG; data_OG{m}.Coefficients.Estimate(c)];
            p{c}(m,1)=ttest_meanSD(data_TM{m}.Coefficients.Estimate(c),data_TM{m}.Coefficients.SE(c),...
                data_OG{m}.Coefficients.Estimate(c),data_OG{m}.Coefficients.SE(c),12);
            
            xlim([.75 2.25])
            ylabel([cond{c}])
            yline(0)
            
            if c==1
                ylim([0 1.5])
            else
                ylim([-.2 1])
            end
            
        end
    end
    
    xlim([.9 2.1])
    xticks([1 2])
    xticklabels(xtl)
    
    
end
[h1,pThreshold1,i11] = BenjaminiHochberg(p{1},0.05);
[h2,pThreshold1,i12] = BenjaminiHochberg(p{2},0.05);


%% Heatmaps 
clear all

groupID='BATR'
adaptation=1; 


if adaptation==1
    load('BAT_adaptation_Indv_0_14-March-2024.mat')
    viewPoints=[440];%PATR - PATS 
else 
    load([groupID, '_post-adaptation_Indv_0_13-March-2024.mat'])
    
end 



muscleOrder={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT',...
    'VM', 'VL', 'RF', 'HIP','TFL', 'GLU'};
ytl=([strcat('f',muscleOrder) strcat('s',muscleOrder)]);  %List of muscle by fast and slow
ytl(end:-1:1) = ytl(:); % We flip the names to match the way labtools grab the data
yt=1:length(ytl); % getting the length of the vector
fs=14;

ex2=[0.2314    0.2980    0.7529];
ex1=[0.7255    0.0863    0.1608];
mid=ones(1,3);
N=100;
gamma=1.5; %gamma > 1 expands the white (mid) part of the map, 'hiding' low values. Gamma<1 does the opposite
gamma=1;
map=[flipud(mid+ (ex1-mid).*([1:N]'/N).^gamma); mid; (mid+ (ex2-mid).*([1:N]'/N).^gamma)];



Y=[];Out=[];
for muscle = 1:28
Y =[Y;model{muscle}.EMGobserved'];
Out =[Out; model{muscle}.EMG_estimated];

end

N=size(Y,2);

% viewPoints=[3,437,443,635]; %PATR - PATS 
binw=4; %Plus minus 2
viewPoints(viewPoints>N-binw/2)=[];
disp(['Data that we are looking at:'])
viewPoints+[-(binw/2):(binw/2)]
Ny=length(viewPoints);
M=length(model);
meanVar=1;%mean(sum(dd.^2,1),2);

figure
for k=1:2
    for i=1:Ny
        switch k
            case 1 % Third row, actual data
                dd=Y(:,viewPoints(i)+[-(binw/2):(binw/2)]);
                nn='data';
            case 2 %Fourth row: one-ahead data predictions
                dd=Out(:,viewPoints(i)+[-(binw/2):(binw/2)]);
                nn={'Data Fit'};
        end
                
                subplot(1,4,k)
                try
                    imagesc(reshape(median(dd,2),12,size(Y,1)/12)')
                    
                catch
                    imagesc(median(dd,2))
                end
                
                set(gca,'XTick',[],'YTick',yt,'YTickLabel',ytl,'FontSize',fs)
                
                ax=gca;
                
                colormap(flipud(map))
                %         caxis([-aC aC])
                caxis([-1 1])
                axis tight
                if k==1
                    title(['Output at t=' num2str(viewPoints(i))])
                    %             txt={'Base Late','Early Adapt','Late Adap','Early Post','Mid Post','Late Post'};
                    %             title(txt{i})
                    ax=gca;
                    ax.Title.FontSize=10;
                end
                if k==2
                    ax=gca;
                    ax.Title.FontSize=10;
                end
                
                
                if i==1
                    ylabel(nn)
                    ax=gca;
                    ax.YAxis.Label.FontWeight='normal';
                    ax.YAxis.Label.FontSize=12;
                end
        end
end
    



for m= 1:length(yt)
    muscle= m; %Needed to skip muscles with VIF>5
    xx=29-m;
    if  VIF_F(muscle,1)<5 % Are early and late R2>0 and VIF<5?
        r2 = my_Rsquared_coeff(median(model{muscle}.EMGobserved(viewPoints(i)+[-(binw/2):(binw/2)],:)), median(model{muscle}.EMG_estimated(:,viewPoints(i)+[-(binw/2):(binw/2)]),2)',1);

        subplot 143
        hold on 
%         plot(m,R2{1}(m,1),'o', "MarkerFaceColor",[0.6350 0.0780 0.1840],"MarkerEdgeColor",[0.6350 0.0780 0.1840])
%         ylabel('R^2')
%         yyaxis left
        scatter(r2,xx,'o','filled',"MarkerFaceColor",[0 0.4470 0.7410])
        xlabel('R^2')

        xline(0)
%         title(xtl{1})
        set(gca,'YTick',yt,'YTickLabel',ytl(end:-1:1),'FontSize',10)
        
        
        subplot 144
        hold on 
%                 yyaxis right
        scatter(VIF_F(m,1),xx,'o',"MarkerEdgeColor",[0.8500 0.3250 0.0980])
        xlabel('VIF')
        set(gca,'YTick',yt,'YTickLabel',ytl(end:-1:1),'FontSize',10)
        xlim([0 3.5])
 
    end
end

 set(gcf,'color','w')