% Grouping muscle by joint bar plots
clear all
load('C3_post-adaptation_04-December-2023.mat')

% 4
[muscle,subj]=find(R2>=0.5); % At the momment of cuttoff values in 50% of variance.
%If the model can reproduce 50% of the variance then it is not a good model

isAffectedSideAll = false(size(R2));
isAffectedSideAll(15:28,1) = true(14,1);
isAffectedSideAll(1:14,2) = true(14,1);
isAffectedSideAll(15:28,3) = true(14,1);
isAffectedSideAll(1:14,4) = true(14,1);
isAffectedSideAll(1:14,5) = true(14,1);
isAffectedSideAll(1:14,6) = true(14,1);
isAffectedSideAll(1:14,7) = true(14,1);
isAffectedSideAll(15:28,8) = true(14,1);
isAffectedSide = false(size(muscle));

for ii = 1:length(muscle)
    isAffectedSide(ii) = isAffectedSideAll(muscle(ii),subj(ii));
end

[muscle_s,subj_s]=find(VIF_F>=5); %finding muscle with colinearity
idx=[];
for i=1:length(muscle_s) %finding the index for muscles with colinearity
    idx=[idx, find(sum(eq([muscle_s(i),subj_s(i)],[muscle,subj]),2)==2)];
end

muscle(idx)=[];subj(idx)=[];  % removing muscles with colinearity
isAffectedSide(idx) = [];

sank_p=[]; sank_n=[];
ship_p=[]; ship_n=[];
sknee_p=[]; sknee_n=[];
fank_p=[]; fank_n=[];
fhip_p=[]; fhip_n=[];
fknee_p=[]; fknee_n=[];

muscleOrder={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT',...
    'VM', 'VL', 'RF', 'HIP','TFL', 'GLU'};
ytl=([strcat('f',muscleOrder) strcat('s',muscleOrder)]);  %List of muscle
ytl(end:-1:1) = ytl(:);
ytl=ytl';
yt=1:length(ytl);

for l=1:length(muscle)
    dataMuscle = [median(reactive_trace(41:45,muscle(l),subj(l)),'omitnan') median(contextual_trace(41:45,muscle(l),subj(l)),'omitnan') muscle(l) subj(l)];
    if muscle(l)<=3 %'sHIP= [reactive contextual muscle_number subject]
        if isAffectedSide(l)
            if isempty(ship_p)
                ship_p(1,:) = dataMuscle;
            else
                ship_p(end+1,:) = dataMuscle;
            end
        else
            if isempty(ship_n)
                ship_n(1,:) = dataMuscle;
            else
                ship_n(end+1,:) = dataMuscle;
            end
        end
    elseif muscle(l)>=4 &&  muscle(l)<=9 %sTHIGH = [reactive contextual muscle_number subject]
        if isAffectedSide(l)
            if isempty(sknee_p)
                sknee_p(1,:) = dataMuscle;
            else
                sknee_p(end+1,:) = dataMuscle;
            end
        else
            if isempty(sknee_n)
                sknee_n(1,:) = dataMuscle;
            else
                sknee_n(end+1,:) = dataMuscle;
            end
        end
    elseif muscle(l)>=10 &&  muscle(l)<=14 %sSHANK = [reactive contextual muscle_number subject]'
        if isAffectedSide(l)
            if isempty(sank_p)
                sank_p(1,:) = dataMuscle;
            else
                sank_p(end+1,:) = dataMuscle;
            end
        else
            if isempty(sank_n)
                sank_n(1,:) = dataMuscle;
            else
                sank_n(end+1,:) = dataMuscle;
            end
        end
    elseif muscle(l)>=15 &&  muscle(l)<=17 %fHIP = [reactive contextual muscle_number subject]'
        if isAffectedSide(l)
            if isempty(fhip_p)
                fhip_p(1,:) = dataMuscle;
            else
                fhip_p(end+1,:) = dataMuscle;
            end
        else
            if isempty(fhip_n)
                fhip_n(1,:) = dataMuscle;
            else
                fhip_n(end+1,:) = dataMuscle;
            end
        end
    elseif muscle(l)>=18 &&  muscle(l)<=23 %fTHIGH = [reactive contextual muscle_number subject]'
        if isAffectedSide(l)
            if isempty(fknee_p)
                fknee_p(1,:) = dataMuscle;
            else
                fknee_p(end+1,:) = dataMuscle;
            end
        else
            if isempty(fknee_n)
                fknee_n(1,:) = dataMuscle;
            else
                fknee_n(end+1,:) = dataMuscle;
            end
        end
    elseif muscle(l)>=24 %fSHANK = [reacrive contextual muscle_number subject]'
        if isAffectedSide(l)
            if isempty(fank_p)
                fank_p(1,:) = dataMuscle;
            else
                fank_p(end+1,:) = dataMuscle;
            end
        else
            if isempty(fank_n)
                fank_n(1,:) = dataMuscle;
            else
                fank_n(end+1,:) = dataMuscle;
            end
        end
    end
end

%% % bar plot by muscle group
close all
colors=["#77AC30";"#7E2F8E";"#D95319"];

for c=1:2
    figure(c);
    subplot 211;
    hold on
    bar([1 4], [median(ship_p(:,c)) median(fhip_p(:,c))],.2,'FaceColor',colors(1))
    errorbar([1 4],[median(ship_p(:,c)) median(fhip_p(:,c))],[iqr(ship_p(:,c)) iqr(fhip_p(:,c))],'.k')
    [p.ship_p(c,1),h.ship_p(c,1),stats.ship_p(c,1)] = signrank(ship_p(:,c));
    [p.fhip_p(c,1),h.fhip_p(c,1),stats.fhip_p(c,1)] = signrank(fhip_p(:,c));
    if p.ship_p(c)<0.05
        plot(1,2,'*r','MarkerSize',15)
    end
    if p.fhip_p(c)<0.05
        plot(4,2,'*r','MarkerSize',15)
    end

    bar([2 5], [median(sknee_p(:,c)) median(fknee_p(:,c))],.2,'FaceColor',colors(2))
    errorbar([2 5],[median(sknee_p(:,c)) median(fknee_p(:,c))],[iqr(sknee_p(:,c)) iqr(fknee_p(:,c))],'.k')
    [p.sknee_p(c,1),h.sknee_p(c,1),stats.sknee_p(c,1)] = signrank(sknee_p(:,c));
    [p.fknee_p(c,1),h.fknee_p(c,1),stats.fknee_p(c,1)] = signrank(fknee_p(:,c));

    if p.sknee_p(c)<0.05
        plot(2,2,'*r','MarkerSize',15);
    end
    if p.fknee_p(c)<0.05
        plot(5,2,'*r','MarkerSize',15)
    end

    bar([3 6], [median(sank_p(:,c)) median(fank_p(:,c))],.2,'FaceColor',colors(3))
    errorbar([3 6],[median(sank_p(:,c)) median(fank_p(:,c))],[iqr(sank_p(:,c)) iqr(fank_p(:,c))],'.k')
    [p.sank_p(c,1),h.sank_p(c,1),stats.sank_p(c,1)] = signrank(sank_p(:,c));
    [p.fank_p(c,1),h.fank_p(c,1),stats.fank_p(c,1)] = signrank(fank_p(:,c));

    if p.sank_p(c)<0.05
        plot(3,2,'*r','MarkerSize',15)
    end
    if p.fank_p(c)<0.05
        plot(6,2,'*r','MarkerSize',15)
    end

    subplot 212;
    hold on
    bar([1 4], [median(ship_n(:,c)) median(fhip_n(:,c))],.2,'FaceColor',colors(1))
    errorbar([1 4],[median(ship_n(:,c)) median(fhip_n(:,c))],[iqr(ship_n(:,c)) iqr(fhip_n(:,c))],'.k')
    [p.ship_n(c,1),h.ship_n(c,1),stats.ship_n(c,1)] = signrank(ship_n(:,c));
    [p.fhip_n(c,1),h.fhip_n(c,1),stats.fhip_n(c,1)] = signrank(fhip_n(:,c));
    if p.ship_n(c)<0.05
        plot(1,2,'*r','MarkerSize',15)
    end
    if p.fhip_n(c)<0.05
        plot(4,2,'*r','MarkerSize',15)
    end

    bar([2 5], [median(sknee_n(:,c)) median(fknee_n(:,c))],.2,'FaceColor',colors(2))
    errorbar([2 5],[median(sknee_n(:,c)) median(fknee_n(:,c))],[iqr(sknee_n(:,c)) iqr(fknee_n(:,c))],'.k')
    [p.sknee_n(c,1),h.sknee_n(c,1),stats.sknee_n(c,1)] = signrank(sknee_n(:,c));
    [p.fknee_n(c,1),h.fknee_n(c,1),stats.fknee_n(c,1)] = signrank(fknee_n(:,c));

    if p.sknee_n(c)<0.05
        plot(2,2,'*r','MarkerSize',15);
    end
    if p.fknee_n(c)<0.05
        plot(5,2,'*r','MarkerSize',15)
    end

    bar([3 6], [median(sank_n(:,c)) median(fank_n(:,c))],.2,'FaceColor',colors(3))
    errorbar([3 6],[median(sank_n(:,c)) median(fank_n(:,c))],[iqr(sank_n(:,c)) iqr(fank_n(:,c))],'.k')
    [p.sank_n(c,1),h.sank_n(c,1),stats.sank_n(c,1)] = signrank(sank_n(:,c));
    [p.fank_n(c,1),h.fank_n(c,1),stats.fank_n(c,1)] = signrank(fank_n(:,c));

    if p.sank_n(c)<0.05
        plot(3,2,'*r','MarkerSize',15)
    end
    if p.fank_n(c)<0.05
        plot(6,2,'*r','MarkerSize',15)
    end
end

color=colormap(turbo(max(subj))) ; % Each participant is a color

for c=1:2
    figure(c);
    subplot 211;
    hold on
    for s=1:max(subj)
        idx_ship_p= find(s==ship_p(:,4));
        idx_fhip_p= find(s==fhip_p(:,4));
        idx_sknee_p= find(s==sknee_p(:,4));
        idx_fknee_p= find(s==fknee_p(:,4));
        idx_sank_p= find(s==sank_p(:,4));
        idx_fank_p= find(s==fank_p(:,4));

        li{s} = scatter([],[],50,'filled','MarkerFaceColor',color(s,:));

        if ~isempty(idx_ship_p)
            li{s}=scatter(1.1,ship_p(idx_ship_p,c),50,'filled','MarkerFaceColor',color(s,:));
        end

        if ~isempty(idx_fhip_p)
            li{s}=scatter(4.1,fhip_p(idx_fhip_p,c),50,'filled','MarkerFaceColor',color(s,:));
        end

        if ~isempty(idx_sknee_p)
            li{s}=scatter(2.1,sknee_p(idx_sknee_p,c),50,'filled','MarkerFaceColor',color(s,:));
        end

        if ~isempty(idx_fknee_p)
            li{s}= scatter(5.1,fknee_p(idx_fknee_p,c),50,'filled','MarkerFaceColor',color(s,:));
        end

        if ~isempty(idx_sank_p)
            li{s}= scatter(3.1,sank_p(idx_sank_p,c),50,'filled','MarkerFaceColor',color(s,:));
        end

        if ~isempty(idx_fank_p)
            li{s}= scatter(6.1,fank_p(idx_fank_p,c),50,'filled','MarkerFaceColor',color(s,:));
        end

        data{s}=li{s}(1);
    end

    xlim([0 7])
    set(gca,'XTick',1:6,'XTickLabel',{'','','','','',''},'FontSize',10)
    xline(3.5)
    ylabel('Affected');

    subplot 212;
    hold on
    for s=1:max(subj)
        idx_ship_n= find(s==ship_n(:,4));
        idx_fhip_n= find(s==fhip_n(:,4));
        idx_sknee_n= find(s==sknee_n(:,4));
        idx_fknee_n= find(s==fknee_n(:,4));
        idx_sank_n= find(s==sank_n(:,4));
        idx_fank_n= find(s==fank_n(:,4));

        li{s} = scatter([],[],50,'filled','MarkerFaceColor',color(s,:));

        if ~isempty(idx_ship_n)
            li{s}=scatter(1.1,ship_n(idx_ship_n,c),50,'filled','MarkerFaceColor',color(s,:));
        end

        if ~isempty(idx_fhip_n)
            li{s}=scatter(4.1,fhip_n(idx_fhip_n,c),50,'filled','MarkerFaceColor',color(s,:));
        end

        if ~isempty(idx_sknee_n)
            li{s}=scatter(2.1,sknee_n(idx_sknee_n,c),50,'filled','MarkerFaceColor',color(s,:));
        end

        if ~isempty(idx_fknee_n)
            li{s}= scatter(5.1,fknee_n(idx_fknee_n,c),50,'filled','MarkerFaceColor',color(s,:));
        end

        if ~isempty(idx_sank_n)
            li{s}= scatter(3.1,sank_n(idx_sank_n,c),50,'filled','MarkerFaceColor',color(s,:));
        end

        if ~isempty(idx_fank_n)
            li{s}= scatter(6.1,fank_n(idx_fank_n,c),50,'filled','MarkerFaceColor',color(s,:));
        end

        data{s}=li{s}(1);
    end

    xlim([0 7])
    set(gca,'XTick',1:6,'XTickLabel',{'sHIP','sKNEE','sANK','fHIP','fKNEE','fANK'},'FontSize',10)
    xline(3.5)
    ylabel('Unaffected');
    legend([data{:}],subID,'Location','best');

    han = axes(gcf,'Visible','off');
    han.XLabel.Visible = 'on';
    han.YLabel.Visible = 'on';
    han.XTickLabel = {'sHIP','sKNEE','sANK','fHIP','fKNEE','fANK'};
    if c==1
        ylabel(han,'W_{reactive}')
    else
        ylabel(han,'W_{contextual}')
    end
    set(gcf,'color','w')
end

