%% Load data and Plot checkerboard for all conditions.
% clear; close all; clc;

% set script parameters, SHOULD CHANGE/CHECK THIS EVERY TIME.
groupID = {'ATS','ATR'};

colorOrder=[[0.4940 0.1840 0.5560];[0.9290 0.6940 0.1250]; [0.4660 0.6740 0.1880]];

scriptDir = fileparts(matlab.desktop.editor.getActiveFilename);


subID = cell(length(groupID), 6);

for g=1:length(groupID)
    files = dir ([scriptDir '/data/' groupID{g} '*params.mat']);
    n_subjects(g) = size(files,1);
    
    
    for i = 1:n_subjects(g)
        sub{g,i} = files(i).name;
        subID{g,i} = sub{g,i}(1:end-10);
        
    end
    
end


subID
 nw=datestr(now,'yy-mm-dd')

%% Norm regressors 

cd(['/Users/dulcemariscal/Documents/GitHub/Generalization_Regressions/RegressionAnalysis/RegModelResults_',nw ,'/IndvResults'])

data=nan(max(n_subjects),length(groupID),6);
for reg=1:2

if reg==1
regressors={'|\DeltaEMG_{adapt}|','|\DeltaEMG_{no-adapt}|','|\DeltaEMG_{env}|'};
tt='Regressors';
c=1:3;
else
regressors={'|\DeltaEMG_{trans1}|','|\DeltaEMG_{trans2}|','|\DeltaEMG_{trans3}|'};
c=4:6;
tt='Transitions';
end



% fh=figure('Units','Normalized','OuterPosition',[0 0 1 1],'NumberTitle', 'off', 'Name',tt);
figure('NumberTitle', 'off', 'Name',tt);


for b=c
    %     subplot(1,6,b)
    hold on
    
    for g=1:length(groupID)
        hold on
        
        if contains(groupID{g},'TR')
            x=[1 4 7 1 4 7];
        else
            x=[2 5 8 2 5 8];
        end
        for s=1:n_subjects(g)
            
            subID{g,s};
            
            load([ subID{g,s} 'defaultsplit_1models_ver00'])
            
            data(s,g,b)= vec_norm(b);
             
            
        end
           
        
        if contains(groupID{g},'TR')
            p{g}=bar(x(b),nanmean(data(:,g,b)),'FaceColor','#0072BD');
       
        else
            p{g}=bar(x(b),nanmean(data(:,g,b)),'FaceColor','w','EdgeColor','#0072BD');

        end
        
        errorbar(x(b),nanmean(data(:,g,b)),nanstd(data(:,g,b))/sqrt(n_subjects(g)),'k')
        
        for  s=1:n_subjects(g)
            
            plot(x(b),data(s,g,b),'*k')
            
        end
    end
    
    
end

legend([p{:}],groupID{:})
ylabel('|\DeltaEMG|')
set(gcf,'color','w');
xticks([1.5 4.5 7.5])
xticklabels(regressors)
end

% box plots 
% figure 
% x1=rmmissing(data(:,1,1));
% x2=rmmissing(data(:,2,1));
% g1 = repmat({'|\DeltaEMG_{adapt1}|'},4,1);
% g2 = repmat({'|\DeltaEMG_{adapt2}|'},7,1);
% x=[x1;x2];
% g=[g1;g2];
% boxplot(x,g,'Notch','on','Whisker',1);


%% Betas bar plots 


cd(['/Users/dulcemariscal/Documents/GitHub/Generalization_Regressions/RegressionAnalysis/RegModelResults_',nw ,'/IndvResults'])
data=[];
for t=1:3
    
    figure
    
    for g=1:length(groupID)
        
         if contains(groupID{g},'TR')
            x=[1 4 7];
        else
            x=[2 5 8];
        end
        
        for b=1:3
            
            for s=1:n_subjects(g)
                
                subID{g,s};
                
                load([ subID{g,s} 'defaultsplit_1models_ver00'])
                
                transition={fitTrans1NoConst, fitTrans2NoConst, fitTrans3NoConst};
                
                hold on
                
                data(s,g,b)=abs(transition{t}.Coefficients.Estimate(b));
                                
            end
            
            if contains(groupID{g},'TR')
                p{b}=bar(x(b),abs(nanmean(data(:,g,b))),'FaceColor',colorOrder(b,:));
            else
                bar(x(b),abs(nanmean(data(:,g,b))),'FaceColor','w','EdgeColor',colorOrder(b,:))
            end
            
            errorbar(x(b),abs(nanmean(data(:,g,b))),abs(nanstd(data(:,g)))/sqrt(n_subjects(g)),'k')
            
            for  s=1:n_subjects(g)
                plot(x(b),data(s,g,b),'*k')
            end
            
        end
    end
    title(['Transition', num2str(t)])
    ylabel('|\beta|')
    xticks([1 2 4 5 7 8])
    xticklabels({'TR','TS','TR','TS','TR','TS'})
    set(gcf,'color','w');
    legend([p{:}],'\beta_{adapt}','\beta_{noadapt}','\beta_{env}')
    %     axis([0 8 -0.1 1.5])
end


%% cd('/Users/dulcemariscal/Documents/GitHub/Generalization_Regressions/RegressionAnalysis/RegModelResults_22-01-31/IndvResults')
data=[];
for t=1:3
    
    figure
    
    for g=1:length(groupID)
        
         if contains(groupID{g},'TR')
            x=[1 4 7];
        else
            x=[2 5 8];
        end
        
        for b=1
            
            for s=1:n_subjects(g)
                
                subID{g,s};
                
                load([ subID{g,s} 'defaultsplit_1models_ver00'])
                
                transition={fitTrans1NoConst, fitTrans2NoConst, fitTrans3NoConst};
                
                hold on
                
                data(s,g,b)=abs(transition{t}.Rsquared.Ordinary);
                                
            end
            
            if contains(groupID{g},'TR')
                p{b}=bar(x(b),abs(nanmean(data(:,g,b))),'FaceColor','#0072BD');
            else
                bar(x(b),abs(nanmean(data(:,g,b))),'FaceColor','w','EdgeColor','#0072BD')
            end
            
            errorbar(x(b),abs(nanmean(data(:,g,b))),abs(nanstd(data(:,g)))/sqrt(n_subjects(g)),'k')
            
            for  s=1:n_subjects(g)
                plot(x(b),data(s,g,b),'*k')
            end
            
        end
    end
    title(['Transition', num2str(t)])
    ylabel('R^2')
    xticks([1 2 4 5 7 8])
    xticklabels({'TR','TS','TR','TS','TR','TS'})
    set(gcf,'color','w');
%     legend([p{:}],'\beta_{adapt}','\beta_{noadapt}','\beta_{env}')
    %     axis([0 8 -0.1 1.5])
end





%% AfterEffects



muscleOrder={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT', 'VM', 'VL', 'RF', 'TFL', 'GLU', 'HIP'};
n_muscles = length(muscleOrder);


ep=defineEpochs_regressionYA('nanmean');
refEp2 = defineReferenceEpoch('OG base',ep);
refEpPost1Early= defineReferenceEpoch('Post1_{Early}',ep);

flip = [1];
norm_ind=nan(max(n_subjects),2);


b=1;
%%
figure
hold on
for g=1:length(groupID)
    
    if contains(groupID{g},'ATR')
        refEp= defineReferenceEpoch('TM base',ep); %fast tied 1 if short split 1, slow tied if 2nd split
        x=[1];
    elseif contains(groupID{g},'ATS')
        refEp= defineReferenceEpoch('OG base',ep);
        x=[2];
    end
    
    
    GroupData=adaptationData.createGroupAdaptData(sub(g,1:n_subjects(g))); %loading the data
    GroupData=GroupData.removeBadStrides; %Removing bad strides
    newLabelPrefix = defineMuscleList(muscleOrder);
    normalizedGroupData = GroupData.normalizeToBaselineEpoch(newLabelPrefix,refEp2); %Normalized by TM base (aka mid baseline)
    
    ll=normalizedGroupData.adaptData{1}.data.getLabelsThatMatch('^Norm');
    l2=regexprep(regexprep(ll,'^Norm',''),'_s','s');
    normalizedGroupData=normalizedGroupData.renameParams(ll,l2);
    newLabelPrefix = regexprep(newLabelPrefix,'_s','s');
    
    
    for s=1:n_subjects(g)
        adaptDataSubject = normalizedGroupData.adaptData{1, s};
        [Data2]=adaptDataSubject.getCheckerboardsData(newLabelPrefix ,refEpPost1Early, refEp,flip);
        norm_ind(s,g)=norm(reshape(Data2,[],1));
        
    end
    
    
    if contains(groupID{g},'TR')
        p{b}=bar(x(b),nanmean(norm_ind(:,g)),'FaceColor','#0072BD');
    else
        bar(x(b),nanmean(norm_ind(:,g)),'FaceColor','w','EdgeColor','#0072BD')
    end
    
    errorbar(x(b),nanmean(norm_ind(:,g)),nanstd(norm_ind(:,g))/sqrt(n_subjects(g)),'k')
    
    for  s=1:n_subjects(g)
        plot(x(b),norm_ind(s,g),'*k')
    end
    
    title(['AfterEffects'])
    ylabel('|\DeltaEMG_{AF}|')
    xticks([1 2])
    xticklabels({'TR','TS'})
    set(gcf,'color','w');
    
end










