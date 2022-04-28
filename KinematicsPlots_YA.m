
% group{1}=adaptationData.createGroupAdaptData({'ATR01params','ATR02params','ATR03params','ATR04params'});
% group{2}=adaptationData.createGroupAdaptData({'ATS02params','ATS03params','ATS04params','ATS05params','ATS06params','ATS07params',...
%      'ATS08params','ATS09params','ATS10params','ATS11params','ATS12params'});
% group{1}=adaptationData.createGroupAdaptData({'ATR01paramsForces','ATR02paramsForces','ATR03paramsForces'});
% group{1}=adaptationData.createGroupAdaptData({'ATS02params','ATS03params','ATS04params','ATS05params','ATS06params'});

% conditions={'OG base','TM slow','TM fast','TM base','Adaptation','Post 1','Post 2','TM mid 1','Pos Short','TM mid 2','Neg Short','TM mid 3'};
% conditions={'OG base','TM slow','TM fast','TM base','Adaptation','Post 1','Post 2','TM mid 1','Pos Short','OG 1','TM mid 3','Neg Short','OG 2'};

% conditions={'OG base','TM slow','TM fast','TM base','Adaptation','Post 1','Post 2','TM mid 1','Pos Short','OG 1','TM mid 2','Neg Short','OG 2'};


% group{1}=adaptationData.createGroupAdaptData({'PATR03params'});
% group{2}=adaptationData.createGroupAdaptData({'PATR03params'});
% conditions={'TM fast','TM slow','TM mid','OG base','TM base','Adaptation','Post 1',...
%     'Post 2','TM mid 1','Pos Short','OG 1','TM mid 2','Neg Short','OG 2',...
%     'Tied','Split 1','Tied 1','Split 2','Tied 3','Split 4',...
%     'Tied 5','Split 6','Tied 7','Split 8','Tied 9','Split 10','Tied 11','Split 12',...
%     'Tied 13','Split 14','Tied 15','Split 16','Tied 17','Split 18','Tied 19','Split 20',...
%     'OG 3','TM mid 3','Neg Short 2','OG 4'};

%  conditions={'TM mid','Pos Short','OG 1','Neg Short','OG 2','TM mid','Adaptation','Post 2',...
%     'Tied','Split 1','Tied 1','Split 2','Tied 3','Split 4',...
%     'Tied 5','Split 6','Tied 7','Split 8','Tied 9','Split 10','Tied 11','Split 12',...
%     'Tied 13','Split 14','Tied 15','Split 16','Tied 17','Split 18','Tied 19','Split 20',...
%     'OG 3'};
% 
% conditions={'Tied Pos','Pos short',...
%     'Tied post Pos','Tied Neg','Neg short','Tied post Neg',...
%     'Tied Ramp','Pos short ramp','Tied post ramp','OG base','TM base',...
%     'Adaptation','Post 1','Post 2','Multiple pos short splits','Split Pos 1',...
%     'Tied Pos 1','Split Pos 2','Tied Pos 3','Split Pos 4','Tied Pos 5','Split Pos 6',...
%     'Tied Pos 7','Split Pos 8','Tied Pos 9','Split Pos 10','Tied Pos 11','Split Pos 12',...
%     'Tied Pos 13','Split Pos 14','Tied Pos 15','Split Pos 16','Tied Pos 17','Split Pos 18',...
%     'Tied Pos 19','Split Pos 20','Tied Pos 21','Split Pos 22','Tied Pos 23','Split Pos 24',...
%     'Tied Pos 25','Split Pos 26','Tied Pos 27','Split Pos 28','Tied Pos 29','Split Pos 30','TM mid 2'};

OA=0;
nimbus=1;
AT=0;
groups=[];
Nimbus=[];
if nimbus==1
    conditions={'OG base','TR base','Adaptation','Post 1'};
    
    load('CTRgroupDataWO_CTR6.mat')
    groups{1}=group;
    control{1}=group;
    load('CTSgroupData.mat')
    groups{2}=group;
    control{2}=group;
    load('NTRgroupData.mat')
    groups{3}=group;
    Nimbus{1}=group;
    load('NTSgroupData.mat')
    groups{4}=group;
    Nimbus{2}=group;
    labels={'CTR','CTS','NTR','NTS'};
    
elseif OA==1
    groups=[];
    load('CgroupData.mat')
    groups{1}=group;
    control{1}=group;
    Nimbus{1}=group;
    
    load('AUFgroupData.mat')
    groups{2}=group;
    control{2}=group;
    
    load('PgroupData.mat')
    groups{3}=group;
    % control{3}=group;
    Nimbus{2}=group;
    labels={'OATR','OATS','Stroke'};
    conditions={'OG base','TM base','Adaptation','Post 1'};
    
elseif AT==1
    
    conditions={'OG base','TM base','Adaptation','Post 1'};
    
    load('ATRgroupData.mat')
    groups{1}=group;
    control{1}=group;
    load('ATSgroupData.mat')
    groups{2}=group;
    control{2}=group;
    labels={'ATR','ATS'};
end
%% 

% conditions={'TM mid','Pos Short','TM mid 2','Neg Short','TM mid 3',...
%     'OG base','TM base','Adaptation','Post 1','Post 2','TM mid 4',...
%     'Multiple Short Splits','Split 1','Tied 1','Split 2','Tied 3',...
%     'Split 4','Tied 5','Split 6','Tied 7','Split 8','Tied 9','Split 10',...
%     'TM mid 5','TM mid 6','Negative adaptation','TM mid 7','Multiple Short  Negative Splits',...
%     'Split Neg 1','Tied Neg 1','Split Neg 2','Tied Neg 3','Split Neg 4','Tied Neg 5','Split Neg 6',...
%     'Tied Neg 7','Split Neg 8','Tied Neg 9','Split Neg 10','TM mid 8'};

% params={'netContributionNorm2'};

params={'UnBiasNormEMG'};
% params={'NormEMG'};
% params={'UnBiasNormEMGasym'};
% params={'NormEMGasym'};
% params={'doubleSupportAsym','doubleSupportDiff','stepTimeContributionNorm2','stanceTimeDiff','stepTimeDiff'};

% params={'stepTimeDiff'};%{'stepTimeContributionNorm2'};%{'doubleSupportAsym'};%,,'doubleSupportDiff','stepTimeDiff'};

poster_colors;
colorOrder=[p_red; p_orange; p_fade_green; p_fade_blue; p_plum; p_green; p_blue; p_fade_red; p_lime; p_yellow; [0 0 0];[0 1 1]];
    
binwidth=5; %Window of the running average
trialMarkerFlag=0; %1 if you want to separete the time course by trial 0 to separece by condition 
indivFlag=0; %0 to plot group mean 1 to plot indiv subjects
indivSubs=[]; %Use when you want to plot a specidfic subject in a group 
% colorOrder=[];%[p_red; p_orange; p_plum;p_fade_green]; %Let the function take care of this at least you wanted in a specific set of color then by my guess and add the list here
biofeedback= 0; % At least that you are providing with biofeedback to the subject
removeBiasFlag=1; %if you want to remove bias 
%%Groups names 
% labels=[];
filterFlag=[0 1]; 
figure 
p=subplot(2,1,1);
plotHandles=p;
alignEnd=0; % # strides align at the end of the trial (PLAY with it as see what happens)
alignIni=0; %  # strides align at the beginning of the trial (PLAY with it as see what happens) 

% adaptData=cellfun(@(x) x.adaptData,groups,'UniformOutput',false); %Notice that adaptDataGroups(1) decide that I only want to plot the CG group 
% [figh,avg,indv]=adaptationData.plotAvgTimeCourse(adaptData,params,conditions,binwidth,trialMarkerFlag,...
%     indivFlag,indivSubs,colorOrder,biofeedback,removeBiasFlag,labels,filterFlag,plotHandles,alignEnd,alignIni);
%%
alignIni=[0 0 120];
alignEnd=[30 40 40];
alignCond={'OGbase','TRbase','Post1'};
axesFontSize=10;
labelFontSize=0;
titleFontSize=12;
row=2;
col=5;

subers=[1:col*row];
subers=(reshape(subers,[col, row]))';

sub=subers(1, 1:4); 


% [ah,figHandle]=optimizedSubPlot(row*col,row,col,'tb',axesFontSize,labelFontSize,titleFontSize);
figure 
p=subplot(2,1,1);
adaptData=cellfun(@(x) x.adaptData,groups,'UniformOutput',false);
plotAvgTimeCourseR01Nimbus(adaptData,params,conditions,binwidth,trialMarkerFlag,indivFlag,...
    indivSubs,colorOrder,0,removeBiasFlag,labels,0,p,alignEnd,alignIni,alignCond, row, col, sub);

legend('AutoUpdate','off')
yline(0)
set(gcf,'color','w');
%%
% figure 
p=subplot(2,3,4);
plotHandles=p;
% if nimbus==1
 conditions={'Post 1'};
alignEnd=40; % # strides align at the end of the trial (PLAY with it as see what happens)
alignIni=120; 
% elseif OA==1
%   conditions={'Washout'};
% end
adaptData=cellfun(@(x) x.adaptData,control,'UniformOutput',false); %Notice that adaptDataGroups(1) decide that I only want to plot the CG group 
[figh,avg,indv]=adaptationData.plotAvgTimeCourse(adaptData,params,conditions,binwidth,trialMarkerFlag,...
    indivFlag,indivSubs,colorOrder,biofeedback,removeBiasFlag,labels,filterFlag,plotHandles,alignEnd,alignIni);
legend('AutoUpdate','off')
yline(0)
% figure 
if nimbus==1 || OA==1
    p=subplot(2,3,5);
    plotHandles=p;
     conditions={'Post 1'};
%     conditions={'Washout'};
if nimbus==1 
    labels2={'NTR','NTS'}; 
    colorOrderNimbus=[p_fade_green; p_fade_blue];
    
elseif OA==1
     labels2={'OATR','Stroke'}; 
     colorOrderNimbus=[p_red; p_fade_green];
end
    
    adaptData=cellfun(@(x) x.adaptData,Nimbus,'UniformOutput',false); %Notice that adaptDataGroups(1) decide that I only want to plot the CG group
    [figh,avg,indv]=adaptationData.plotAvgTimeCourse(adaptData,params,conditions,binwidth,trialMarkerFlag,...
        indivFlag,indivSubs,colorOrderNimbus,biofeedback,removeBiasFlag,labels2,filterFlag,plotHandles,alignEnd,alignIni);
    legend('AutoUpdate','off')
    yline(0)
end
%%


if nimbus==1
%     epochNames={'OG_{l}','TR_{l}','Adapt_{e}','Adapt_{l}','P1_{e}','P1_{l}'};
%     conditions={'OG base','TR base','Adaptation','Adaptation','Post 1','Post 1'}; %Change conditions names to your own!
    epochNames={'P1_{e}','P1_{l}'};
    conditions={'Post 1','Post 1'};
    
else
    epochNames={'OG_{l}','TM_{l}','Adapt_{e}','Adapt_{l}','P1_{e}','P1_{l}'};
    conditions={'OG base','TM base','Adaptation','Adaptation','Post 1','Post 1'}; %Change conditions names to your own!
end
% strideNo=[-20 -20 5 -20 5 -20 ]; %Positive vaues define inital; negative values define # strides at end of that condition
strideNo=[5 -40 ]; 
% exemptFirst=[ 1 1 1 1 1 1 ]; %Number of strides you want to ignore at the beginning of the condition
exemptFirst=[ 1 1 ]; 
% exemptLast=[ 5 5 0 5 0 5]; %Number of strides you want to ignore at the end of the condition
exemptLast=[0 5];
summaryMethod={'nanmean'}; %Method to analyze bar plots 
summaryMethod={'nanmedian'}; %Method to analyze bar plots 
shortName=[];

[epochs] = defineEpochs(epochNames,conditions,strideNo,exemptFirst,exemptLast,summaryMethod,shortName);

% [data,validStrides,everyStrideData]=getEpochData(groups{1},epochs,params,0); %This function output a matrix with the data from the epochs FYI

% params={'stepLengthAsym'}; %you can plot which ever parameter you are interested on check into the adaptData.data for labels
% binwidth=5; %Window of the running average
% trialMarkerFlag=0; %1 if you want to separete the time course by trial 0 to separece by condition 
indivFlag=1; %0 to plot group mean 1 to plot indiv subjects
% indivSubs=[]; %Use when you want to plot a specidfic subject in a group 
% % colorOrder=[]; %Let the function take care of this at least you wanted in a specific set of color then by my guess and add the list here
% biofeedback= 0; % At least that you are providing with biofeedback to the subject
% % labels={'CTR','CTS','NTR','NTS'}; %Groups names 
% % removeBaseEpochFlag=1;
% removeBiasFlag=1;
% alignEnd=0;
significanceThreshold=[];
posthocGroupFlag=[];
posthocEpochFlag=[];
% figure
% p=subplot(1,1,1);
p=subplot(2,3,6);
plotHandles=p;
posthocGroupByEpochFlag=[];
posthocEpochByGroupFlag=[];
medianFlag=[];

[allData]=groupAdaptationData.plotMultipleEpochBars(groups,params,epochs,indivFlag,labels,plotHandles,colorOrder,...
    medianFlag,significanceThreshold,posthocGroupFlag,posthocEpochFlag,posthocGroupByEpochFlag,...
    posthocEpochByGroupFlag,removeBiasFlag);
% [fh,ph,allData]=adaptationData.plotGroupedTimeAndEpochBars(groups,params,epochs,binwidth,trialMarkerFlag,indivFlag,indivSubs,colorOrder,biofeedback,labels,medianFlag,removeBaseEpochFlag,alignEnd,significanceThreshold,posthocGroupFlag,posthocEpochFlag,posthocGroupByEpochFlag,posthocEpochByGroupFlag);
%%

OA=0;
nimbus=1;
AT=0;
groups=[];
Nimbus=[];
if nimbus==1
 conditions={'OG base','TR base','Adaptation','Post 1'};

load('CTRgroupDataWO_CTR6.mat')
groups{1}=group;
control{1}=group;
% load('CTSgroupData.mat')
% groups{2}=group;
% control{2}=group;
% load('NTRgroupData.mat')
% groups{3}=group;
% Nimbus{1}=group;
% load('NTSgroupData.mat')
% groups{1}=group;
% Nimbus{2}=group;
 labels={'CTR','CTS','NTR','NTS'}; 
elseif OA==1
groups=[];
load('CgroupData.mat')
groups{1}=group;
control{1}=group;
Nimbus{1}=group;

load('AUFgroupData.mat')
groups{2}=group;
control{2}=group;

load('PgroupData.mat')
groups{3}=group;
% control{3}=group;
Nimbus{2}=group;
labels={'OATR','OATS','Stroke'}; 


% 
 conditions={'OG base','TM base','Adaptation','Post 1'};
 
elseif AT==1
    
     conditions={'OG base','TM base','Adaptation','Post 1','Post 2'};

load('ATRgroupData.mat')
groups{1}=group;
control{1}=group;
% load('ATSgroupData.mat')
% groups{2}=group;
% control{2}=group;
%    labels={'ATR','ATS'};  
end
% conditions={'TM mid','Pos Short','TM mid 2','Neg Short','TM mid 3',...
%     'OG base','TM base','Adaptation','Post 1','Post 2','TM mid 4',...
%     'Multiple Short Splits','Split 1','Tied 1','Split 2','Tied 3',...
%     'Split 4','Tied 5','Split 6','Tied 7','Split 8','Tied 9','Split 10',...
%     'TM mid 5','TM mid 6','Negative adaptation','TM mid 7','Multiple Short  Negative Splits',...
%     'Split Neg 1','Tied Neg 1','Split Neg 2','Tied Neg 3','Split Neg 4','Tied Neg 5','Split Neg 6',...
%     'Tied Neg 7','Split Neg 8','Tied Neg 9','Split Neg 10','TM mid 8'};

% params={'netContributionNorm2'};

params={'UnBiasNormEMG'};
% params={'NormEMG'};
% params={'UnBiasNormEMGasym'};
% params={'NormEMGasym'};
% params={'doubleSupportAsym','doubleSupportDiff','stepTimeContributionNorm2','stanceTimeDiff','stepTimeDiff'};

% params={'stepTimeDiff'};%{'stepTimeContributionNorm2'};%{'doubleSupportAsym'};%,,'doubleSupportDiff','stepTimeDiff'};

poster_colors;
colorOrder=[p_red; p_orange; p_fade_green; p_fade_blue; p_plum; p_green; p_blue; p_fade_red; p_lime; p_yellow; [0 0 0];[0 1 1]];
    
binwidth=5; %Window of the running average
trialMarkerFlag=0; %1 if you want to separete the time course by trial 0 to separece by condition 
indivFlag=1; %0 to plot group mean 1 to plot indiv subjects
indivSubs=[]; %Use when you want to plot a specidfic subject in a group 
% colorOrder=[];%[p_red; p_orange; p_plum;p_fade_green]; %Let the function take care of this at least you wanted in a specific set of color then by my guess and add the list here
biofeedback= 0; % At least that you are providing with biofeedback to the subject
removeBiasFlag=1; %if you want to remove bias 
%%Groups names 
% labels=[];
filterFlag=[]; 
figure 
p=subplot(1,1,1);
plotHandles=p;
alignEnd=0; % # strides align at the end of the trial (PLAY with it as see what happens)
alignIni=0; %  # strides align at the beginning of the trial (PLAY with it as see what happens) 

adaptData=cellfun(@(x) x.adaptData,groups,'UniformOutput',false); %Notice that adaptDataGroups(1) decide that I only want to plot the CG group 
[figh,avg,indv]=adaptationData.plotAvgTimeCourse(adaptData,params,conditions,binwidth,trialMarkerFlag,...
    indivFlag,indivSubs,colorOrder,biofeedback,removeBiasFlag,labels,filterFlag,plotHandles,alignEnd,alignIni);
legend('AutoUpdate','off')
yline(0)
set(gcf,'color','w');
%%
% groups{1}=adaptationData.createGroupAdaptData({'ST10params','ST11params','ST13params','ST14params','ST16params','ST19params'});
% 
% 
% conditions={'TM base','Adaptation','TM post'};
% params={'stanceTimeDiff'};%{'stepTimeContributionNorm2'};%{'doubleSupportAsym'};%,,,'stepTimeDiff'};
% 
% % params={'doubleSupportAsym','doubleSupportDiff','stepTimeContributionNorm2','stanceTimeDiff','stepTimeDiff'};
% poster_colors;
% colorOrder=[p_red; p_orange; p_fade_green; p_fade_blue; p_plum; p_green; p_blue; p_fade_red; p_lime; p_yellow; [0 0 0];[0 1 1]];
%     
% binwidth=5; %Window of the running average
% trialMarkerFlag=0; %1 if you want to separete the time course by trial 0 to separece by condition 
% indivFlag=0; %0 to plot group mean 1 to plot indiv subjects
% indivSubs=[]; %Use when you want to plot a specidfic subject in a group 
% colorOrder=[];%[p_red; p_orange; p_plum;p_fade_green]; %Let the function take care of this at least you wanted in a specific set of color then by my guess and add the list here
% biofeedback= 0; % At least that you are providing with biofeedback to the subject
% removeBiasFlag=1; %if you want to remove bias 
% labels=[]; %Groups names 
% filterFlag=[]; 
% plotHandles=[];
% alignEnd=0; % # strides align at the end of the trial (PLAY with it as see what happens)
% alignIni=0; %  # strides align at the beginning of the trial (PLAY with it as see what happens) 
% 
% adaptData=cellfun(@(x) x.adaptData,groups,'UniformOutput',false); %Notice that adaptDataGroups(1) decide that I only want to plot the CG group 
% [figh,avg,indv]=adaptationData.plotAvgTimeCourse(adaptData,params,conditions,binwidth,trialMarkerFlag,...
%     indivFlag,indivSubs,colorOrder,biofeedback,removeBiasFlag,labels,filterFlag,plotHandles,alignEnd,alignIni);
% set(gcf,'color','w'); 