clear
clc
close all

%% Definitions
% cd('C:\Users\gonza\OneDrive - University of Pittsburgh\Desktop\Data');
cd('C:\Users\mag356\Desktop\R01Nimbus2021\Data');

param={'netContributionNorm2'};
trialMarkerFlag=[0 0 0];
IndivSubList=[];
conds1={'OG Base', 'TR Base', 'Adaptation', 'Post 1'};
% conds1={'OG Base','TR Base','Post 1'};
epochs={'EarlyPost1', 'LatePost1'};
poster_colors;
colorOrder=[];
indivSubFlag=[]; %Change this only if you have multiple groups
plotOnly1Sub='';
removeBias=1;

useFirstTied=false;

binWidth=5;

alignIni=150;
alignEnd=40;
alignCond={'Post1'};

% alignIni=[];
% alignEnd=[];
% alignCond=[];

%% General definitions (do not change)

if useFirstTied
    CTS=adaptationData.createGroupAdaptData({'CTS_03_TRbaseFirst','CTS_04_TRbaseFirst','CTS_05_TRbaseFirst','CTS_06_TRbaseFirst','CTS_07_TRbaseFirst'});
    CTR=adaptationData.createGroupAdaptData({'CTR_02_TRbaseFirst','CTR_03_TRbaseFirst','CTR_04_TRbaseFirst','CTR_05_TRbaseFirst','CTR_06_TRbaseFirst'});
else
    CTS=adaptationData.createGroupAdaptData({'CTS_03','CTS_04','CTS_05','CTS_06','CTS_07'});
    CTR=adaptationData.createGroupAdaptData({'CTR_02','CTR_03','CTR_04','CTR_05','CTR_06'});
end
NTR=adaptationData.createGroupAdaptData({'NTR_01','NTR_02','NTR_03','NTR_04','NTR_05'});
NTS=adaptationData.createGroupAdaptData({'NTS_01','NTS_02','NTS_04','NTS_06','NTS_07'}); % 'NTS_08' device malfunction,  'NTS_03','NTS_05' and 'NTS_09' removed because least amount of AE in Post1

StudyData{1}=CTS;
StudyData{2}=CTR;
StudyData{3}=NTS;
StudyData{4}=NTR;

StudyData=cellfun(@(x) x.adaptData,StudyData,'UniformOutput',false);
groups={'Control Testing', 'Control Training', 'Nimbus Testing', 'Nimbus Training'};
%% Set up the physical Figure

axesFontSize=10;
labelFontSize=0;
titleFontSize=12;
row=2;
col=5;

[ah,figHandle]=optimizedSubPlot(row*col,row,col,'tb',axesFontSize,labelFontSize,titleFontSize);


%% Timecourse plots Control Testing

subers=[1:col*row];
subers=(reshape(subers,[col, row]))';

sub=subers(1, 1:4); 

plotAvgTimeCourseR01Nimbus(StudyData,param,conds1,binWidth,trialMarkerFlag,indivSubFlag,IndivSubList,colorOrder,0,removeBias,groups,0,ah,alignEnd,alignIni,alignCond, row, col, sub);
hold on;
xlimits=xlim;
xlimits=[xlimits(1):1:xlimits(2)];
p3=plot(xlimits, zeros(1,length(xlimits)), 'k', 'LineWidth', 1.5);
% legend off;
hold off;

if removeBias
    title('Step Length Asymmetry Unbiased');
else
    title('Step Length Asymmetry Biased');
end
ylim([-0.2 0.2]);


%% Zoom in on the Time courses

sub=subers(2, 1:2); 

plotAvgTimeCourseR01Nimbus({StudyData{3:4}},param,{'Post 1'},binWidth,trialMarkerFlag,indivSubFlag,IndivSubList,[p_fade_green; p_fade_blue],0,removeBias,{'Nimbus Testing', 'Nimbus Training'},0,ah,alignEnd,alignIni,alignCond, row, col, sub);
hold on;
xlimits=xlim;
xlimits=[xlimits(1):1:xlimits(2)];
p3=plot(xlimits, zeros(1,length(xlimits)), 'k', 'LineWidth', 1.5);
hold off;

% ylim([-0.025 0.2]);


sub=subers(2, 3:4); 

plotAvgTimeCourseR01Nimbus({StudyData{1:2}},param,{'Post 1'},binWidth,trialMarkerFlag,indivSubFlag,IndivSubList,colorOrder,0,removeBias,{'Control Testing', 'Control Training'},0,ah,alignEnd,alignIni,alignCond, row, col, sub);
hold on;
xlimits=xlim;
xlimits=[xlimits(1):1:xlimits(2)];
p3=plot(xlimits, zeros(1,length(xlimits)), 'k', 'LineWidth', 1.5);
hold off;

% ylim([-0.025 0.2]);


%% 

R=getResults_R01NimV2(StudyData,param,groups,removeBias); % Unbiased


%Epoch 1
barGroupsSingleR01NimbusV2(StudyData,R,groups,param,{'EarlyPost1'},colorOrder, row, col,[5],1);
% ylim([-0.025 0.2]);

barGroupsSingleR01NimbusV2(StudyData,R,groups,param,{'LatePost1'},colorOrder, row, col,[10],1);
% ylim([-0.25 0.2]);



