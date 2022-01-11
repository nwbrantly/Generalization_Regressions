
group{1}=adaptationData.createGroupAdaptData({'ATR01params','ATR02params','ATR03params','ATR04params'});
% group{1}=adaptationData.createGroupAdaptData({'ATR01paramsForces','ATR02paramsForces','ATR03paramsForces'});
group{1}=adaptationData.createGroupAdaptData({'ATS02params','ATS03params','ATS04params','ATS05params',});

conditions={'OG base','TM slow','TM fast','TM base','Adaptation','Post 1','Post 2','TM mid 1','Pos Short','TM mid 2','Neg Short','TM mid 3'};

conditions={'OG base','TM slow','TM fast','TM base','Adaptation','Post 1','Post 2','TM mid 1','Pos Short','OG 1','TM mid 3','Neg Short','OG 2'};
% 
% conditions={'TM base','Adaptation','Post 1'};

params={'netContributionNorm2'};


% params={'doubleSupportAsym','doubleSupportDiff','stepTimeContributionNorm2','stanceTimeDiff','stepTimeDiff'};

% params={'stepTimeDiff'};%{'stepTimeContributionNorm2'};%{'doubleSupportAsym'};%,,'doubleSupportDiff','stepTimeDiff'};

poster_colors;
colorOrder=[p_red; p_orange; p_fade_green; p_fade_blue; p_plum; p_green; p_blue; p_fade_red; p_lime; p_yellow; [0 0 0];[0 1 1]];
    
binwidth=5; %Window of the running average
trialMarkerFlag=0; %1 if you want to separete the time course by trial 0 to separece by condition 
indivFlag=0; %0 to plot group mean 1 to plot indiv subjects
indivSubs=[]; %Use when you want to plot a specidfic subject in a group 
colorOrder=[];%[p_red; p_orange; p_plum;p_fade_green]; %Let the function take care of this at least you wanted in a specific set of color then by my guess and add the list here
biofeedback= 0; % At least that you are providing with biofeedback to the subject
removeBiasFlag=0; %if you want to remove bias 
labels=[]; %Groups names 
filterFlag=[]; 
plotHandles=[];
alignEnd=0; % # strides align at the end of the trial (PLAY with it as see what happens)
alignIni=0; %  # strides align at the beginning of the trial (PLAY with it as see what happens) 

adaptData=cellfun(@(x) x.adaptData,group,'UniformOutput',false); %Notice that adaptDataGroups(1) decide that I only want to plot the CG group 
[figh,avg,indv]=adaptationData.plotAvgTimeCourse(adaptData,params,conditions,binwidth,trialMarkerFlag,...
    indivFlag,indivSubs,colorOrder,biofeedback,removeBiasFlag,labels,filterFlag,plotHandles,alignEnd,alignIni);
set(gcf,'color','w');

%%


group{1}=adaptationData.createGroupAdaptData({'ST10params','ST11params','ST13params','ST14params','ST16params','ST19params'});


conditions={'TM base','Adaptation','TM post'};
params={'stanceTimeDiff'};%{'stepTimeContributionNorm2'};%{'doubleSupportAsym'};%,,,'stepTimeDiff'};

% params={'doubleSupportAsym','doubleSupportDiff','stepTimeContributionNorm2','stanceTimeDiff','stepTimeDiff'};
poster_colors;
colorOrder=[p_red; p_orange; p_fade_green; p_fade_blue; p_plum; p_green; p_blue; p_fade_red; p_lime; p_yellow; [0 0 0];[0 1 1]];
    
binwidth=5; %Window of the running average
trialMarkerFlag=0; %1 if you want to separete the time course by trial 0 to separece by condition 
indivFlag=0; %0 to plot group mean 1 to plot indiv subjects
indivSubs=[]; %Use when you want to plot a specidfic subject in a group 
colorOrder=[];%[p_red; p_orange; p_plum;p_fade_green]; %Let the function take care of this at least you wanted in a specific set of color then by my guess and add the list here
biofeedback= 0; % At least that you are providing with biofeedback to the subject
removeBiasFlag=1; %if you want to remove bias 
labels=[]; %Groups names 
filterFlag=[]; 
plotHandles=[];
alignEnd=0; % # strides align at the end of the trial (PLAY with it as see what happens)
alignIni=0; %  # strides align at the beginning of the trial (PLAY with it as see what happens) 

adaptData=cellfun(@(x) x.adaptData,group,'UniformOutput',false); %Notice that adaptDataGroups(1) decide that I only want to plot the CG group 
[figh,avg,indv]=adaptationData.plotAvgTimeCourse(adaptData,params,conditions,binwidth,trialMarkerFlag,...
    indivFlag,indivSubs,colorOrder,biofeedback,removeBiasFlag,labels,filterFlag,plotHandles,alignEnd,alignIni);
set(gcf,'color','w');