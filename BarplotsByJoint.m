% Grouping muscle by joint bar plots
clear all
% load('BATS_post-adaptationAdaptation_01-December-2023.mat')
% load('CTS_post-adaptation_04-December-2023.mat')
% load('BATR_post-adaptation_29-November-2023.mat')
load('AUF_post-adaptation_30-November-2023.mat')

% 4
[muscle,subj]=find(R2>=0.5); % At the momment of cuttoff values in 50% of variance.
%If the model can reproduce 50% of the variance then it is not a godo model

[muscle_s,subj_s]=find(VIF_F>=5); %finding muscle with colinearity
idx=[];
for i=1:length(muscle_s) %finding the index for muscles with colinearity

    idx=[idx, find(sum(eq([muscle_s(i),subj_s(i)],[muscle,subj]),2)==2)];
end

muscle(idx)=[];subj(idx)=[];  % removing muscles with colinearity

sank=[];
ship=[];
sknee=[];
fank=[];
fhip=[];
fknee=[];

muscleOrder={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT',...
    'VM', 'VL', 'RF', 'HIP','TFL', 'GLU'};
ytl=([strcat('f',muscleOrder) strcat('s',muscleOrder)]);  %List of muscle
ytl(end:-1:1) = ytl(:);
ytl=ytl';
yt=1:length(ytl);

for l=1:length(muscle)

    if muscle(l)<=3 %'sHIP= [reacrive contextual muscle_number subject]

        if isempty(ship)==1
            ship(1,:)=[median(reactive_trace(41:45,muscle(l),subj(l)),'omitnan') median(contextual_trace(41:45,muscle(l),subj(l)),'omitnan') muscle(l) subj(l)];

        else
            ship(end+1,:)=[ median(reactive_trace(41:45,muscle(l),subj(l)),'omitnan') median(contextual_trace(41:45,muscle(l),subj(l)),'omitnan') muscle(l) subj(l)];
        end

    elseif muscle(l)>=4 &&  muscle(l)<=9 %sTHIGH = [reacrive contextual muscle_number subject]

        if isempty(sknee)==1
            sknee(1,:)=[median(reactive_trace(41:45,muscle(l),subj(l)),'omitnan') median(contextual_trace(41:45,muscle(l),subj(l)),'omitnan') muscle(l) subj(l)];
        else
            sknee(end+1,:)=[median(reactive_trace(41:45,muscle(l),subj(l)),'omitnan') median(contextual_trace(41:45,muscle(l),subj(l)),'omitnan') muscle(l) subj(l)];
        end


    elseif muscle(l)>=10 &&  muscle(l)<=14 %sSHANK = [reacrive contextual muscle_number subject]'

        if isempty(sank)==1
            sank(1,:)=[median(reactive_trace(41:45,muscle(l),subj(l)),'omitnan') median(contextual_trace(41:45,muscle(l),subj(l)),'omitnan') muscle(l) subj(l)];
        else
            sank(end+1,:)=[median(reactive_trace(41:45,muscle(l),subj(l)),'omitnan') median(contextual_trace(41:45,muscle(l),subj(l)),'omitnan') muscle(l) subj(l)];
        end


    elseif muscle(l)>=15 &&  muscle(l)<=17 %fHIP = [reacrive contextual muscle_number subject]'
        if isempty(fhip)==1

            fhip(1,:)=[median(reactive_trace(41:45,muscle(l),subj(l)),'omitnan') median(contextual_trace(41:45,muscle(l),subj(l)),'omitnan') muscle(l) subj(l)];
        else
            fhip(end+1,:)=[ median(reactive_trace(41:45,muscle(l),subj(l)),'omitnan') median(contextual_trace(41:45,muscle(l),subj(l)),'omitnan') muscle(l) subj(l)];
        end


    elseif muscle(l)>=18 &&  muscle(l)<=23 %fTHIGH = [reacrive contextual muscle_number subject]'
        if isempty(fknee)==1
            fknee(1,:)=[median(reactive_trace(41:45,muscle(l),subj(l)),'omitnan') median(contextual_trace(41:45,muscle(l),subj(l)),'omitnan') muscle(l) subj(l)];
        else
            fknee(end+1,:)=[median(reactive_trace(41:45,muscle(l),subj(l)),'omitnan') median(contextual_trace(41:45,muscle(l),subj(l)),'omitnan') muscle(l) subj(l)];
        end

    elseif muscle(l)>=24 %fSHANK = [reacrive contextual muscle_number subject]'
        if isempty(fank)==1
            fank(1,:)=[median(reactive_trace(41:45,muscle(l),subj(l)),'omitnan') median(contextual_trace(41:45,muscle(l),subj(l)),'omitnan') muscle(l) subj(l)];
        else
            fank(end+1,:)=[median(reactive_trace(41:45,muscle(l),subj(l)),'omitnan') median(contextual_trace(41:45,muscle(l),subj(l)),'omitnan') muscle(l) subj(l)];
        end

    end
end



%% % bar plot by muscle group
close all
colors=["#77AC30";"#7E2F8E";"#D95319"];


for c=1:2
    figure(c)
    hold on
    bar([1 4], [median(ship(:,c)) median(fhip(:,c))],.2,'FaceColor',colors(1))
    errorbar([1 4],[median(ship(:,c)) median(fhip(:,c))],[iqr(ship(:,c)) iqr(fhip(:,c))] ,'.k')
    [p.ship(c,1),h.ship(c,1),stats.ship(c,1)] = signrank(ship(:,c));
    [p.fhip(c,1),h.fhip(c,1),stats.fhip(c,1)] = signrank(fhip(:,c));
    if p.ship(c)<0.05
        plot(1,2,'*r','MarkerSize',15)
        %         text(1,2.05,num2str(round(p.ship(c),2)))
    end

    if p.fhip(c)<0.05
        plot(4,2,'*r','MarkerSize',15)
        %         text(4,2.05,num2str(round(p.fhip(c),2)))
    end
    %     plot(1.1,ship(:,c),'*k')
    %     plot(4.1,fhip(:,c),'*k')

    bar([2 5], [median(sknee(:,c)) median(fknee(:,c))],.2,'FaceColor',colors(2))
    errorbar([2 5],[median(sknee(:,c)) median(fknee(:,c))],[iqr(sknee(:,c)) iqr(fknee(:,c))] ,'.k')
    [p.sknee(c,1),h.sknee(c,1),stats.sknee(c,1)] = signrank(sknee(:,c));
    [p.fknee(c,1),h.fknee(c,1),stats.fknee(c,1)] = signrank(fknee(:,c));

    if p.sknee(c)<0.05
        plot(2,2,'*r','MarkerSize',15);
        %         text(2,2.05,num2str(round(p.sknee(c),2)))

    end

    if p.fknee(c)<0.05
        plot(5,2,'*r','MarkerSize',15)
        %         text(5,2.05,num2str(round(p.fknee(c),2)))
    end

    %     plot(2.1,sknee(:,c),'*k')
    %     plot(5.1,fknee(:,c),'*k')

    hold on
    bar([3 6], [median(sank(:,c)) median(fank(:,c))],.2,'FaceColor',colors(3))
    errorbar([3 6],[median(sank(:,c)) median(fank(:,c))],[iqr(sank(:,c))/sqrt(length(sank)) iqr(fank(:,c))/sqrt(length(sank)) ] ,'.k')
    [p.sank(c,1),h.sank(c,1),stats.sank(c,1)] = signrank(sank(:,c));
    [p.fank(c,1),h.fank(c,1),stats.fank(c,1)] = signrank(fank(:,c));

    if p.sank(c)<0.05
        plot(3,2,'*r','MarkerSize',15)
        %         text(3,2.05,num2str(round(p.sank(c),2)))
    end

    if p.fank(c)<0.05
        plot(6,2,'*r','MarkerSize',15)
        %         text(6,2.05,num2str(round(p.fank(c),2)))
    end

    %     plot(3.1,sank(:,c),'*k')
    %     plot(6.1,fank(:,c),'*k')

    xlim([0 7])
    set(gca,'XTick',[1:6],'XTickLabel',{'sHIP','sKNEE','sANK','fHIP','fKNEE','fANK'},'FontSize',10)
    xline(3.5)

    if c==1
        ylabel('W_{reactive}')
    else
        ylabel('W_{contextual}')

    end
    set(gcf,'color','w')
end

color=colormap(turbo(max(subj))) ; % Each participant is a color

for c=1:2
    figure(c);
    hold on
    for s=1:max(subj)
        idx_ship= find(s==ship(:,4));
        idx_fhip= find(s==fhip(:,4));
        idx_sknee= find(s==sknee(:,4));
        idx_fknee= find(s==fknee(:,4));
        idx_sank= find(s==sank(:,4));
        idx_fank= find(s==fank(:,4));

        li{s} = scatter([],[],50,'filled','MarkerFaceColor',color(s,:));

        if ~isempty(idx_ship)
            li{s}=scatter(1.1,ship(idx_ship,c),50,'filled','MarkerFaceColor',color(s,:));
            %             plot(1.1,ship(idx_ship,c),'.','MarkerFaceColor',color(s,:),'MarkerSize',25)
        end

        if ~isempty(idx_fhip)
            li{s}=scatter(4.1,fhip(idx_fhip,c),50,'filled','MarkerFaceColor',color(s,:));
            %             plot(4.1,fhip(idx_fhip,c),'.','MarkerFaceColor',color(s,:),'MarkerSize',25)
        end

        if ~isempty(idx_sknee)
            li{s}=scatter(2.1,sknee(idx_sknee,c),50,'filled','MarkerFaceColor',color(s,:));
            %             plot(1.1,ship(idx_ship,c),'.','MarkerFaceColor',color(s,:),'MarkerSize',25)
        end

        if ~isempty(idx_fknee)
            li{s}= scatter(5.1,fknee(idx_fknee,c),50,'filled','MarkerFaceColor',color(s,:));
            %             plot(4.1,fhip(idx_fhip,c),'.','MarkerFaceColor',color(s,:),'MarkerSize',25)
        end

        if ~isempty(idx_sank)
            li{s}= scatter(3.1,sank(idx_sank,c),50,'filled','MarkerFaceColor',color(s,:));
            %             plot(1.1,ship(idx_ship,c),'.','MarkerFaceColor',color(s,:),'MarkerSize',25)
        end

        if ~isempty(idx_fank)
            li{s}= scatter(6.1,fank(idx_fank,c),50,'filled','MarkerFaceColor',color(s,:));
            %             plot(4.1,fhip(idx_fhip,c),'.','MarkerFaceColor',color(s,:),'MarkerSize',25)
        end

        data{s}=li{s}(1);
    end
    legend([data{:}],subID,'Location','best');
end


%% Separating by anterior and posterior muscles 

data_slow= [ship; sknee; sank];
data_fast= [fhip; fknee; fank];

% anterior = {'TA', 'VM', 'VL', 'RF', 'TFL','HIP'};
% posterior = {'LG','MG','SOL','PER','SEMT','SEMB','BF', 'GLU'};

% groupong by flexors and extensor
extensors={'LG','MG','SOL','PER','VM','VL','RF',};
flexors = {'TA','SEMT','SEMB','BF','HIP','TFL','GLU'};

anterior = extensors;
posterior= flexors;

anterior_list=([strcat('f',anterior) strcat('s',anterior)]);  %List of muscle
anterior_list(end:-1:1) =anterior_list(:);
anterior_list=anterior_list';

posterior_list =([strcat('f',posterior) strcat('s',posterior)]);  %List of muscle
posterior_list(end:-1:1) =posterior_list(:);
posterior_list=posterior_list';

anterior_muscles_slow=[];
anterior_muscles_fast=[];
for ant= 1:length(anterior_list)
    
idx= find(strcmp( [ytl(:)], anterior_list{ant} ));

idx_muscle_s= find(idx==data_slow(:,3));
idx_muscle_f= find(idx==data_fast(:,3));

anterior_muscles_slow= [anterior_muscles_slow;data_slow(idx_muscle_s,:) ];
anterior_muscles_fast= [anterior_muscles_fast;data_fast(idx_muscle_f,:) ];
end

posterior_muscles_slow=[];
posterior_muscles_fast=[];
for ant= 1:length(posterior_list)
    
idx= find(strcmp( [ytl(:)], posterior_list{ant} ));

idx_muscle_s= find(idx==data_slow(:,3));
idx_muscle_f= find(idx==data_fast(:,3));

posterior_muscles_slow= [posterior_muscles_slow;data_slow(idx_muscle_s,:) ];
posterior_muscles_fast= [posterior_muscles_fast;data_fast(idx_muscle_f,:) ];

end

colors=["#77AC30";"#7E2F8E";"#D95319"];
for c=1:2
    figure(c)
    hold on
    bar([1 3], [median(anterior_muscles_slow(:,c)) median(anterior_muscles_fast(:,c))],.2,'FaceColor',colors(1))
    errorbar([1 3],[median(anterior_muscles_slow(:,c)) median(anterior_muscles_fast(:,c))],[iqr(anterior_muscles_slow(:,c)) iqr(anterior_muscles_fast(:,c))] ,'.k')
    [p.anteriorslow(c,1),h.anteriorslow(c,1),] = signrank(anterior_muscles_slow(:,c));
    [p.anteriorfast(c,1),h.fhip(c,1),stats.anteriorfast(c,1)] = signrank(anterior_muscles_fast(:,c));
    if p.anteriorslow(c)<0.05
        plot(1,2,'*r','MarkerSize',15)
        %         text(1,2.05,num2str(round(p.ship(c),2)))
    end

    if p.anteriorfast(c)<0.05
        plot(3,2,'*r','MarkerSize',15)
        %         text(4,2.05,num2str(round(p.fhip(c),2)))
    end
    %     plot(1.1,ship(:,c),'*k')
    %     plot(4.1,fhip(:,c),'*k')

    bar([2 4], [median(posterior_muscles_slow(:,c)) median(posterior_muscles_fast(:,c))],.2,'FaceColor',colors(2))
    errorbar([2 4],[median(posterior_muscles_slow(:,c)) median(posterior_muscles_fast(:,c))],[iqr(posterior_muscles_slow(:,c)) iqr(posterior_muscles_fast(:,c))] ,'.k')
    [p.posterior_muscles_slow(c,1),h.posterior_muscles_slow(c,1),stats.posterior_muscles_slow(c,1)] = signrank(posterior_muscles_slow(:,c));
    [p.posterior_muscles_fast(c,1),h.posterior_muscles_fast(c,1),stats.posterior_muscles_fast(c,1)] = signrank(posterior_muscles_fast(:,c));

    if p.posterior_muscles_slow(c)<0.05
        plot(2,2,'*r','MarkerSize',15);
        %         text(2,2.05,num2str(round(p.sknee(c),2)))

    end

    if p.posterior_muscles_fast(c)<0.05
        plot(4,2,'*r','MarkerSize',15)
        %         text(5,2.05,num2str(round(p.posterior_muscles_fast(c),2)))
    end


    xlim([0 5])
%     set(gca,'XTick',[1:4],'XTickLabel',{'slow_{anterior}','fast_{anterior}','slow_{posterior}','fast_{posterior}'},'FontSize',10)
  set(gca,'XTick',[1:4],'XTickLabel',{'slow_{extensor}','fast_{extensor}','slow_{flexor}','fast_{flexor}'},'FontSize',10)
xline(2.5)

    if c==1
        ylabel('W_{reactive}')
    else
        ylabel('W_{contextual}')

    end
    set(gcf,'color','w')
end


color=colormap(turbo(max(subj))) ; % Each participant is a color

for c=1:2
    figure(c);
    hold on
    for s=1:max(subj)
        idx_anterior_muscles_slow= find(s==anterior_muscles_slow(:,4));
        idx_anterior_muscles_fast= find(s==anterior_muscles_fast(:,4));
        idx_posterior_muscles_slow= find(s==posterior_muscles_slow(:,4));
        idx_posterior_muscles_fast= find(s==posterior_muscles_fast(:,4));


        li{s} = scatter([],[],50,'filled','MarkerFaceColor',color(s,:));

        if ~isempty(idx_anterior_muscles_slow)
            li{s}=scatter(1.1,anterior_muscles_slow(idx_anterior_muscles_slow,c),50,'filled','MarkerFaceColor',color(s,:));
            %             plot(1.1,anterior_muscles_slow(idx_anterior_muscles_slow,c),'.','MarkerFaceColor',color(s,:),'MarkerSize',25)
        end

        if ~isempty(idx_anterior_muscles_fast)
            li{s}=scatter(3.1,anterior_muscles_fast(idx_anterior_muscles_fast,c),50,'filled','MarkerFaceColor',color(s,:));
            %             plot(4.1,anterior_muscles_fast(idx_anterior_muscles_fast,c),'.','MarkerFaceColor',color(s,:),'MarkerSize',25)
        end

        if ~isempty(idx_posterior_muscles_slow)
            li{s}=scatter(2.1,posterior_muscles_slow(idx_posterior_muscles_slow,c),50,'filled','MarkerFaceColor',color(s,:));
            %             plot(1.1,anterior_muscles_slow(idx_anterior_muscles_slow,c),'.','MarkerFaceColor',color(s,:),'MarkerSize',25)
        end

        if ~isempty(idx_posterior_muscles_fast)
            li{s}= scatter(4.1,posterior_muscles_fast(idx_posterior_muscles_fast,c),50,'filled','MarkerFaceColor',color(s,:));
            %             plot(4.1,anterior_muscles_fast(idx_anterior_muscles_fast,c),'.','MarkerFaceColor',color(s,:),'MarkerSize',25)
        end


        data{s}=li{s}(1);
    end
    legend([data{:}],subID,'Location','best');
end
