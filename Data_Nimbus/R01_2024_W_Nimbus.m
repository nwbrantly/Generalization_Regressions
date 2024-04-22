% Plotting the value of the initial aftereffect for the reactive and
% contextual process. We are eliminating muscles with VIF> 5 (colinearily)
% and R2<0.2 (bad reconstruction)  


steps=41:45;
groupID={'CTR','NTR','VATR','CTS','NTS','VATS'};

figure
hold on
reactive=[];
contextual=[];

colors= [[0.8500 0.3250 0.0980];[0 0.4470 0.7410];[0.4660 0.6740 0.1880];[1 1 1];[1 1 1];[1 1 1]];
colorEdge= [[1 1 1];[1 1 1];[1 1 1];[0.8500 0.3250 0.0980];[0 0.4470 0.7410];[0.4660 0.6740 0.1880]];

for g=1:length(groupID) %loop across groups
    
    load([groupID{g},'_post-adaptation_22-February-2024.mat'])
    reconstruction_contextual=squeeze(mean(contextual_trace(steps,:,:),'omitnan')); %Early aftereffect for the contextual 
    reconstruction_reactive=squeeze(mean(reactive_trace(steps,:,:),'omitnan')); %Early aftereffect for the reactive
    
    for s=1:size(reconstruction_contextual,2) %Loop across participants marking those muscles with VIF>5 and R2<0.2 as NaN
        reconstruction_reactive(find(VIF_F(:,s)>5),s)=nan; 
        reconstruction_reactive(find(R2{s}(:,1)<0.2),s)=nan;
        reconstruction_contextual(find(VIF_F(:,s)>5),s)=nan;
        reconstruction_contextual(find(R2{s}(:,1)<0.2),s)=nan;
        reconstruction_reactive(find(reconstruction_reactive(:,s)<0),s)=nan;
        reconstruction_contextual(find(reconstruction_contextual(:,s)<0),s)=nan;
        
    end
     %reorgnanize the data as a vectors 
    reactive(:,g)= reshape(reconstruction_reactive,28*5,1);
    contextual(:,g)= reshape(reconstruction_contextual,28*5,1); 
    
    % Addign the number of muscle that pass our criteria for reactive &
    % contextual 
    n_muscle(g)=sum(~isnan(reactive(:,g)));
    n_muscle_C(g)=sum(~isnan(contextual(:,g)));
    
    
    %Plot reactive data 
    subplot 211;hold on
    scatter(g,reactive(:,g),'o','filled','MarkerFaceColor',[.7 .7 .7]) %Plotting the individual muscle that pass the criteria  
    %Plotting mean and standard desviation 
    plot(g+.1,mean(reactive(:,g),'omitnan'),'sk','MarkerFaceColor',colors(g,:),...
        'MarkerSize',15,'MarkerEdgeColor',colorEdge(g,:))
    errorbar(g+.1,mean(reactive(:,g),'omitnan'),std(reactive(:,g),'omitnan'),'k')
    
    h= ttest(reactive(:,g)); %testing if the distribution is different from zero 
    
    
    if h==1 % If the distribution is difference from zero. We plot a magenta * on top 
        plot(g, 2, '*m') 
    end
    
    %Plot contextual data 
    subplot 212 ;hold on
    scatter(g,contextual(:,g),'o','filled','MarkerFaceColor',[.7 .7 .7]) %Plotting the individual muscle that pass the criteria  
    plot(g+.1,mean(contextual(:,g),'omitnan'),'sk','MarkerFaceColor',colors(g,:),...
        'MarkerSize',15,'MarkerEdgeColor',colorEdge(g,:))
    errorbar(g+.1,mean(contextual(:,g),'omitnan'),std(contextual(:,g),'omitnan'),'k')
    h= ttest(contextual(:,g));
    
    if h==1 % If the distribution is difference from zero. We plot a magenta * on top 
        plot(g, 2, '*m')
    end
    
    
end

% Update x-label names to match the groups
subplot 211;hold on
set(gca,'XTick',1:6,'XTickLabel',groupID,'FontSize',10)
yline(0)
ylabel({'Reactive';'(A.U)'})

subplot 212 ;hold on
set(gca,'XTick',1:6,'XTickLabel',groupID,'FontSize',10)
yline(0)
ylabel({'Contextual';'(A.U)'})
set(gcf,'color','w')
%% Number of muscle per group that pass our criteria VIF>5 and R2<0.2 - Needs input from aboce section 
figure
hold on

for g=1:6
    
    bar(g,n_muscle(g),'FaceColor',colors(g,:),'EdgeColor',colorEdge(g,:))
    
end
title('Total of muscle per group')
ylabel('Number of muscles')
set(gca,'XTick',1:6,'XTickLabel',groupID,'FontSize',10)




%% Number of muscle with reactive and contextual > 0 - Needs inputs from the section 1 
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



%% ttest between groups, you need to pic the groups - Needs inputs from the section 1

[h_reactive,p_reactive] = ttest2(reactive(:,4),reactive(:,5))
[h_contextual,p_contextual] = ttest2(contextual(:,4),contextual(:,5))

%% Two-way ANOVA - Needs inputs from the section 1

 %Grouping data by adaptation condition 
CT=[contextual(:,1);contextual(:,4)];
NT=[contextual(:,2);contextual(:,5)];
VAT=[contextual(:,3);contextual(:,6)];

g1 = [ones(1,28*5) 2*ones(1,28*5) 3*ones(1,28*5) 1*ones(1,28*5) 2*ones(1,28*5) 3*ones(1,28*5) ];
g2=[ones(1,28*5*3) 2*ones(1,28*5*3)];

% reactive
c=reshape(reactive,28*5*6,1) ; %reshaping the data as a vector 
[~,~,stats] = anovan(c,{g1 g2},"Model","interaction", ... %Using anovan since we have nan on the data 
    "Varnames",["Adaptation","PostAdaptation"]);

%contextual
c=reshape(contextual,28*5*6,1) ; %reshaping the data as a vector 
[~,~,stats] = anovan(c,{g1 g2},"Model","interaction", ... %Using anovan since we have nan on the data 
    "Varnames",["Adaptation","PostAdaptation"]);
%% Stroke Data 
%% Data from C3S01 with and without the nimbus 
%Same plot as above the costum part is the groups colors 

steps=41:45; %Steps that we want to look at 
groupID={'C3TR_sub1','MWSTR_sub1','C3TS_sub1','MWSTS_sub1'}; %Grpups name 
figure
colors= [[1 1 1];[1 1 1];[0.8500 0.3250 0.0980];[0 0.4470 0.7410]]; %color per group 
colorEdge= [[0.8500 0.3250 0.0980];;[0 0.4470 0.7410];[1 1 1];[1 1 1];]; %Colors per group 


figure
hold on
reactive=[];
contextual=[];

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

%% First 9 participants of the C3 study data 

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

%% Number of muscles that are difference from zero 
%- Using the SE form the regression, if cross zero than we do not count that muscle 
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
%%
clear all
binwith=5
analysis=0;isF=0;

% Define which group you want to plot
% groupID={'VATR','VATS'};
% groupID={'NTR','NTS'};
 groupID={'CTR','CTS'};
% groupID={'BATR','BATS'};
f=[1:14 1:14;1:14 1:14];

for g=1:2 %Group loop (Device or W/O device)
    
    %Load data
    if strncmpi(groupID{g},'BAT',3) %See if the first letters of the group are BAT
        load([groupID{g},'_post-adaptation_Indv_0_13-March-2024.mat'])
    else
        load([groupID{g},'_post-adaptation_Indv_0_08-March-2024.mat'])
    end
    
    
    for muscle=1:length(labels) %Define muscles to plot (This data set has a total of 28 muscles)
        
        r2 = my_Rsquared_coeff(model{muscle}.EMGobserved', model{muscle}.EMG_estimated,1); % Compute the R2 throuhgtout the data
        figure(f(g,muscle))
        
        % Pick the data that you want to plot
        Wdata(:,1)=[reactive_trace(41:end,muscle)]'; %Data organization and selecting only the postadaptation data
        Wdata(:,2)=[contextual_trace(41:end,muscle)]'; %Data organization and selecting only the postadaptation data
        r2=[r2(:,41:end)];
        
        % Making indeces where r2 is negative nana
        idx= find(r2<0); % finding when the r2 is negative
        Wdata(idx,1)=nan; %reactive
        Wdata(idx,2)=nan; %contextual
        r2(:,idx)=nan;
        
        % Getting the average of the data while ignoting nan
        avg(:,1)=movmean(Wdata(:,1),binwith,'omitnan'); %reactive
        avg(:,2)=movmean(Wdata(:,2),binwith,'omitnan'); %contextual
        r2_avg= movmean(r2,binwith,'omitnan');
        
        % Alternative: we get the average and then we mark data as NaN.
        % This leads to more empty data points
        %             idx= find(r2_avg<0);
        %             avg(idx,:)=nan;
        %             r2_avg(:,idx)=nan;
        
        
        %Plot reactive time course 
        if muscle<15 %Choosing where to plot the time course 
            subplot(size(Wdata,2)+1,2,1)
        else
            subplot(size(Wdata,2)+1,2,2)
        end
        
        hold on
        scatter(1:length(avg),avg(:,1),'filled','DisplayName',groupID{g}); %plot  
        legend %Adding the groupID 
        title(labels(muscle).Data) %adding muscle that we are looking at on the title 
        yline(0,'HandleVisibility','off')
        ylabel({'Reactive';'(A.U)'})
        xlabel('strides')
        
        
        %Plot  contextual time course 
        if size(Wdata,2)>=2
            if muscle<15
                subplot(size(Wdata,2)+1,2,3)
            else
                subplot(size(Wdata,2)+1,2,4)
            end
            hold on
            scatter(1:length(avg),avg(:,2),'filled','DisplayName',groupID{g}); % plot 
            legend  %Adding the groupID 
            yline(0,'HandleVisibility','off')
            ylabel({'Contextual';'(A.U)'})
            xlabel('strides')
        end
        

        
        % Plot the time course of the R2
        if muscle<15
            subplot(size(Wdata,2)+1,2,5)
        else
            subplot(size(Wdata,2)+1,2,6)
        end
        
        hold on
        plot(1:length(r2_avg),r2_avg,'DisplayName',groupID{g}); %plot 
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


%% Time courses and R2 - Organizaton similar to Adwoa's paper
%%

binwith=5;
analysis=0;isF=0;

%Organizing the data 
groupID{1}={'VATR','VATS'}; % Very low contextual similarity 
groupID{2}={'BATR','BATS'}; % Low contextual similarity 
groupID{3}={'NTR','NTS'}; %High contextual similarity 


colors= [[1 1 1];[0.4660 0.6740 0.1880]; ...
    [1 1 1]; [0 0.4470 0.7410];...
    [1 1 1];[0.8500 0.3250 0.0980]]; %Defining colors for the marker 
colorEdge=[[0.4660 0.6740 0.1880];[0 0 0];...
    [0 0.4470 0.7410];[0 0 0];...
    [0.8500 0.3250 0.0980];[0 0 0]]; %Defining colors for the marker edge 
i=0;

 f=[1 5 9]; % hard code position of the plots 

for muscle=[20]
    counter=0; % This counter is use to determine the location of the plot 
    i=0; % this counter is use for calling the color to be use 
   
    figure %Everytime that we plot a new muscle we create a figure 
    
    
    for level=1:3
       
        counter=counter+1;   % This counter is use to determine the location of the plot 
         
        for g=1:2 % Groups (W device and W/O device) 
            i=i+1; 
            
            if level==2 % loading the data - I keep track of the date to make sure I know wish changes were made 
                load([groupID{level}{g},'_post-adaptation_Indv_0_13-March-2024.mat'])
            else
                load([groupID{level}{g},'_post-adaptation_Indv_0_08-March-2024.mat'])
            end
            labels(muscle).Data %Display for the use to know which muscle we are plotting
            groupID{level}{g} %Display the name of the group for the user to know where we are in the loop 
            
            
            r2 = my_Rsquared_coeff(model{muscle}.EMGobserved', model{muscle}.EMG_estimated,1); %Calculating the R2 time course 
            
            % Pick the data that you want to plot
            Wdata(:,1)=[reactive_trace(44:150,muscle)]'; %Data organization and selecting only the postadaptation data 
            Wdata(:,2)=[contextual_trace(44:150,muscle)]'; %Data organization and selecting only the postadaptation data 
            r2=[r2(:,44:150)];

            % Making indeces where r2 is negative nana 
            idx= find(r2<0); % finding when the r2 is negative 
            Wdata(idx,1)=nan; %reactive
            Wdata(idx,2)=nan; %contextual 
            r2(:,idx)=nan;
            model{muscle}.EMG_estimated(:,idx)=nan;
            model{muscle}.EMGobserved(idx,:)=nan;
            
            % Getting the average of the data while ignoting nan 
            avg(:,1)=movmean(Wdata(:,1),binwith,'omitnan'); %reactive
            avg(:,2)=movmean(Wdata(:,2),binwith,'omitnan'); %contextual 
            r2_avg= movmean(r2,binwith,'omitnan');
            
            % Alternative: we get the average and then we mark data as NaN.
            % This leads to more empty data points
%             idx= find(r2_avg<0);
%             avg(idx,:)=nan;
%             r2_avg(:,idx)=nan;
            
            %getting confidance interval 
            temp=  movmean(model{muscle}.EMGobserved(44:150,:)',binwith,2,'omitnan');
            aftereffects=temp(:,1) ;
            temp=  movmean(model{muscle}.EMG_estimated(:,44:150),binwith,2,'omitnan');
            estimated =temp(:,1) ;
            R2 = my_Rsquared_coeff(aftereffects,estimated,1);
            m= fitlm(model{muscle}.C,aftereffects,'VarNames',{'Reactive','Contextual', labels(muscle).Data(1:end-1)},'Intercept',false);
             
            
            % Ploting the reactive time course 
            subplot(3,4,f(counter))          
            hold on
            scatter(1:length(avg),avg(:,1),'filled','DisplayName',groupID{level}{g},'MarkerEdgeColor',colorEdge(i,:),'MarkerFaceColor',colors(i,:)); %"#77AC30" )%
            legend
            title(labels(muscle).Data)            
            yline(0,'HandleVisibility','off')
            ylabel({'Reactive';'(A.U)'});xlabel('strides')
            if muscle==6
                ylim([-.2 .7])
            elseif muscle==20
                ylim([-.5 .8])
            end
            
            
            % Plotting the reactive box data. Note we are plottign the mean
            % of the first 5 strides 
            subplot(1,4,2);hold on
            errorbar(g+counter/10,avg(1,1),m.Coefficients.SE(1),'k')
            plot(g+counter/10,avg(1,1),'s','MarkerEdgeColor',colorEdge(i,:),'MarkerFaceColor',colors(i,:),'MarkerSize',15)
            
            xlim([0.8 2.5])
            ylabel({'Reactive';'(A.U)'})
            

            
            if size(Wdata,2)>=2
                % Ploting the contextual time course 
                subplot(3,4,f(counter)+2)               
                hold on
                scatter(1:length(avg), avg(:,2),'filled','DisplayName',groupID{level}{g},'MarkerEdgeColor',colorEdge(i,:),'MarkerFaceColor',colors(i,:));
                legend
                yline(0,'HandleVisibility','off')
                ylabel({'Contextual';'(A.U)'})
                xlabel('strides')
                title(labels(muscle).Data)
                if muscle==6
                    ylim([-.2 .7])
                elseif muscle==20
                    ylim([-1 .8])
                    
                end
            % Plotting the contextual box data. Note we are plottign the mean
            % of the first 5 strides                 
                subplot(1,4,4);hold on
                errorbar(g+counter/10,avg(1,2),m.Coefficients.SE(2),'k')
                plot(g+counter/10,avg(1,2),'s','MarkerEdgeColor',colorEdge(i,:),'MarkerFaceColor',colors(i,:),'MarkerSize',15)
                xlim([0.8 2.5])
                ylabel({'Contextual';'(A.U)'})
                
            end
               
        end
        
    end
    set(gcf,'color','w') %setting figure background as white 
    
    % Hard code x-labels and y-limits for the box plot figures 
    subplot(1,4,2);hold on
    set(gca,'XTick',1:2,'XTickLabel',{'W Device','W/O Device'},'FontSize',10)
    title('Early Post-adaptation')
    if muscle==6
        ylim([-.2 .7])
    elseif muscle==20
        ylim([-.2 .8])
    end
    
    
    subplot(1,4,4);hold on
    set(gca,'XTick',1:2,'XTickLabel',{'W Device','W/O Device'},'FontSize',10)
    title('Early Post-adaptation')
    if muscle==6
        ylim([-.2 .8])
    elseif muscle==20
        ylim([-1 .8])
        
    end
    
    set(gcf,'renderer','painters') % Rendering data to make it editable in illustraitor 
end


%% EMG TRACES - DATA FROMG YOUNG ADULTS BATR03

% Set muscle to plot
% muscle={ 'TFL', 'GLU','HIP', 'SEMB', 'SEMT','BF', 'VM', 'VL', 'RF','SOL', 'LG', 'MG','TA', 'PER'}; %muscles that you want to plot 

muscle={ 'TFL'}; %muscles that you want to plot 
normalize = 1;  % 1 to normalize data
normCond = {'TM mid 1'}; % Condition that you want to use to normalize the data 

conds={'TM mid 1','Pos Short'};
late=[1 0];  %0 average of initial strides 1 average of the last strides
strides=[10 30]; % Number of strides that you are going to average 
IgnoreStridesEarly=[1 30]; %number of strides that you are going to ignore at the beginning
plotEMGtraces(expData,conds,muscle,late,strides,normalize,normCond,IgnoreStridesEarly); 


conds={'TM mid 1','Multiple Pos Shorts Splits'};
late=[1 1];  %0 average of initial strides 1 average of the last strides
strides=[10 40]; % Number of strides that you are going to average 
IgnoreStridesEarly=[1 0]; %number of strides that you are going to ignore at the beginning
plotEMGtraces(expData,conds,muscle,late,strides,normalize,normCond,IgnoreStridesEarly); 


conds={'TM mid 1','Adaptation'};
late=[1 1];  %0 average of initial strides 1 average of the last strides
strides=[10 40]; % Number of strides that you are going to average 
IgnoreStridesEarly=[1 0]; %number of strides that you are going to ignore at the beginning
plotEMGtraces(expData,conds,muscle,late,strides,normalize,normCond,IgnoreStridesEarly);







