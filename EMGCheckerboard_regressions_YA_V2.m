%This file will run EMG regressions and plot checkerboards.
% Focus on subjects that performed 2 sets of
% pos/neg short perturbations
%
% 1)load data; assuming the data is saved under currentDir/data/
% 2)plot checkboards for all relevant epochs, save figures, can be turned
% off by setting plotAllEpoch to false
% 3)plot checkboards for regression related epoch (regressors), save figures,
%run regression and save the model results.
% - can plot and run regression for both indidual subjects or group subjects (set by plotGroup flag or when there are more than 1 subjects provided),
% turn off individual subjects plotting by setting to false
% The results are saved under
% currentDir/RegressionAnalysis/RegModelResults_V##. If there are code
% changes that's worth a version update, search for _V## and then update
% the version number to avoid overwrite.

%% Load data and Plot checkerboard for all conditions.
clear; close all; clc;

% set script parameters, SHOULD CHANGE/CHECK THIS EVERY TIME.
groupID = 'ATS';
saveResAndFigure = true;
plotAllEpoch = true;
plotIndSubjects = true;
plotGroup = true;
bootstrap=true;

scriptDir = fileparts(matlab.desktop.editor.getActiveFilename);
files = dir ([scriptDir '/data/' groupID '*params.mat']);

n_subjects = size(files,1);
subID = cell(1, n_subjects);
sub=cell(1,n_subjects);
for i = 1:n_subjects
    sub{i} = files(i).name;
    subID{i} = sub{i}(1:end-10);
   
end
subID

regModelVersion =  'default';
%% load and prep data
GroupData=adaptationData.createGroupAdaptData(sub); %loading the data
GroupData=GroupData.removeBadStrides; %Removing bad strides

muscleOrder={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT', 'VM', 'VL', 'RF', 'TFL', 'GLU', 'HIP'};

n_muscles = length(muscleOrder);
n_subjects = length(subID);

ep=defineEpochs_regressionYA('nanmean');
refEp = defineReferenceEpoch('TM base',ep);
refEpLate = defineReferenceEpoch('Adaptation',ep);
% refEpSlow = defineReferenceEpoch('TM slow',ep);
refEp2 = defineReferenceEpoch('OG base',ep);

newLabelPrefix = defineMuscleList(muscleOrder);
normalizedGroupData = GroupData.normalizeToBaselineEpoch(newLabelPrefix,refEp2); %Normalized by TM base (aka mid baseline)

ll=normalizedGroupData.adaptData{1}.data.getLabelsThatMatch('^Norm');


l2=regexprep(regexprep(ll,'^Norm',''),'_s','');
normalizedGroupData=normalizedGroupData.renameParams(ll,l2);

newLabelPrefix = regexprep(newLabelPrefix,'_s','');

%% Removal of muscle with porblems during data collection
% Do not remove mouscle for group analysis
if contains(groupID,'ATR')
    badMuscleNames = {'sGLUs','fGLUs',};
elseif contains(groupID,'ATS02')
    badMuscleNames = {'sRFs','fRFs','sHIPs','fHIPs'};
elseif contains(groupID,'ATS')
    badMuscleNames = {'sRFs','fRFs'};
end
badMuscleIdx=[];
for bm = badMuscleNames
    badMuscleIdx = [badMuscleIdx, find(ismember(newLabelPrefix,bm))];
end
newLabelPrefix = newLabelPrefix(setdiff(1:end, badMuscleIdx))
%% Plot epochs
ep=defineEpochs_regressionYA('nanmean');
plotGroup=1;
plotIndSubjects = 0;
if plotGroup
    summaryflag='nanmedian';
        adaptDataSubject = normalizedGroupData;
    figureSaveId = groupID;
    sub=1;
    
elseif plotIndSubjects
    summaryflag=[];
        adaptDataSubject = normalizedGroupData.adaptData{1, i};
    figureSaveId = subID{i};
    sub=n_subjects;
end

if plotAllEpoch
    for i = 1:sub
                
        fh=figure('Units','Normalized','OuterPosition',[0 0 1 1]);
        ph=tight_subplot(1,length(ep)-1,[.03 .005],.04,.04);
        flip=true;
        
        refEp=[];
        adaptDataSubject.plotCheckerboards(newLabelPrefix,ep(1:8,:),fh,ph(1,1:8),refEp,flip,summaryflag);%plot all epochs normalized by the mid baseline
        adaptDataSubject.plotCheckerboards(newLabelPrefix,ep(10:17,:),fh,ph(1,9:length(ep)-1),refEp,flip,summaryflag);%plot all epochs normalized by the mid baseline
       
        
        set(ph(:,1),'CLim',[-1 1]*1);
        set(ph(:,2:end),'YTickLabels',{},'CLim',[-1 1]*1);
        set(ph,'FontSize',8)
        pos=get(ph(1,end),'Position');
        axes(ph(1,end))
        colorbar
        set(ph(1,end),'Position',pos);
        set(gcf,'color','w');
        
        if (saveResAndFigure)
            resDir = [scriptDir '/RegressionAnalysis/RegModelResults_V14/'];
            if not(isfolder(resDir))
                mkdir(resDir)
            end
            saveas(fh, [resDir subID{i} '_AllEpochCheckerBoard.png'])
            saveas(fh, [resDir subID{i} '_AllEpochCheckerBoard'],'epsc')
        end
    end
end

%% Prepare data for regressor checkerboards and regression model


usefft = 0; %Turn on to use the transpose of the EMG(+) as the estimated for EMG(-)
normalizeData = 0; %Vector lenght normalization

%type of regression that we want to run
regModelVersion='default';

% In case that we have multiple cases. Let see if this needs to be updated
splitCount=1;

%For group data try to use the median to be consistnace with Pablo's

summaryflag='nanmean';

% Defining the epochs of interest
ep=defineEpochs_regressionYA(summaryflag);

%Getting the names for specific epocjhs
OGbase=defineReferenceEpoch('OG base',ep);
TMbase= defineReferenceEpoch('TM base',ep);
TMbeforePos= defineReferenceEpoch('TM_{beforePos}',ep);
PosShort= defineReferenceEpoch('PosShort_{early}',ep);
TMbeforeNeg=defineReferenceEpoch('TM_{beforeNeg}',ep);
NegShort=defineReferenceEpoch('NegShort_{early}',ep);
Pos1_Early=defineReferenceEpoch('Post1_{Early}',ep);
AdaptLate=defineReferenceEpoch('Adaptation',ep);
Pos1_Late=defineReferenceEpoch('Post1_{Late}',ep);
Pos2_Early=defineReferenceEpoch('Post2_{Early}',ep);

if (contains(groupID, 'ATR'))
    AfterPos= defineReferenceEpoch('TM_{afterPos}',ep);
    AfterNeg=defineReferenceEpoch('TM_{afterNeg}',ep);
    
elseif (contains(groupID, 'ATS'))
    AfterPos= defineReferenceEpoch('OG_{afterPos}',ep);
    AfterNeg=defineReferenceEpoch('OG_{afterNeg}',ep);
end


%% plot checkerboard and run regression per subject
if plotIndSubjects
    for i = 1:n_subjects
        
        for flip=1:2
            
            [Data,regressorNames,fh]=RegressionsGeneralization(newLabelPrefix,normalizedGroupData,[],1,0,NegShort,TMbeforeNeg,PosShort,TMbeforePos,...
                AdaptLate,Pos1_Early,Pos1_Late,Pos2_Early, AfterPos, OGbase, TMbase, i,flip);
            
            nw=datestr(now,'yy-mm-dd');
            
            resDir = [scriptDir '/RegressionAnalysis/RegModelResults_', nw,'/IndvResults/'];
            if (saveResAndFigure)
                if not(isfolder(resDir))
                    mkdir(resDir)
                end
                saveas(fh, [resDir subID{i} '_Checkerboard_ver' num2str(usefft) num2str(normalizeData) regModelVersion 'split_' num2str(splitCount) '.png'])
                saveas(fh, [resDir subID{i} '_Checkerboard_ver' num2str(usefft) num2str(normalizeData) regModelVersion 'split_' num2str(splitCount)],'epsc')
            end
            
            if flip ~= 2
                % run regression and save results
                format compact % format loose %(default)
                % not normalized first, then normalized, arugmnets order: (Data, normalizeData, isGroupData, dataId, resDir, saveResAndFigure, version, usefft)
                fprintf('\n')
                splitCount
                runRegression_Generalization(Data, false, false, [subID{i} regModelVersion 'split_' num2str(splitCount)], resDir, saveResAndFigure, regModelVersion, usefft)
                runRegression_Generalization(Data, true, false, [subID{i} regModelVersion 'split_' num2str(splitCount)], resDir, saveResAndFigure, regModelVersion, usefft)
                
                asymCos = findCosBtwAsymOfEpochs(Data, size(newLabelPrefix,2))
            end
        end
        
    end
end
%% plot checkerboards and run regression per group
if length(subID) > 1 || plotGroup
    
    for  flip = 2
        [Data,regressorNames]=RegressionsGeneralization(newLabelPrefix,normalizedGroupData,[],0,1,NegShort,TMbeforeNeg,PosShort,TMbeforePos,...
            AdaptLate,Pos1_Early,Pos1_Late,Pos2_Early, AfterPos, OGbase, TMbase,[],flip);
        
        
        nw=datestr(now,'yy-mm-dd');
        resDir = [scriptDir '/RegressionAnalysis/RegModelResults_',nw ,'/GroupResults/'];
        if (saveResAndFigure)
            if not(isfolder(resDir))
                mkdir(resDir)
            end
%             saveas(fh, [resDir groupID '_group_Checkerboard_ver' num2str(usefft) num2str(normalizeData) regModelVersion 'split_' num2str(splitCount) 'flip_' num2str(flip) '.png'])
%             saveas(fh, [resDir groupID '_group_Checkerboard_ver' num2str(usefft) num2str(normalizeData) regModelVersion 'split_' num2str(splitCount) 'flip_' num2str(flip)], 'epsc')
        end
        
        if flip ~= 2 %run regression on the full (not asymmetry) data
            % run regression and save results
            format compact % format loose %(default)
            % not normalized first, then normalized, arugmnets order: (Data, normalizeData, isGroupData, dataId, resDir, saveResAndFigure, version, usefft) 
            runRegression_Generalization(Data, false, true, [groupID regModelVersion 'split_' num2str(splitCount) 'flip_' num2str(flip)], resDir, saveResAndFigure, regModelVersion, usefft)
            runRegression_Generalization(Data, true, true, [groupID regModelVersion 'split_' num2str(splitCount) 'flip_' num2str(flip)], resDir, saveResAndFigure, regModelVersion, usefft)
        else
            asymCos = findCosBtwAsymOfEpochs(Data, size(newLabelPrefix,2))
        end
    end
end
    
%% EMG aftereffects

refEpPost1Early= defineReferenceEpoch('Post1_{Early}',ep);
if contains(groupID,'ATR')
    refEp= defineReferenceEpoch('TM base',ep); %fast tied 1 if short split 1, slow tied if 2nd split
elseif contains(groupID,'ATS')
    refEp= defineReferenceEpoch('OG base',ep);
end

fh=figure('Units','Normalized','OuterPosition',[0 0 1 1]);
ph=tight_subplot(1,n_subjects+1,[.03 .005],.04,.04);
flip = [1];

nw=datestr(now,'yy-mm-dd');
% nw='22-01-20';
if plotIndSubjects
    for i = 1:n_subjects
        adaptDataSubject = normalizedGroupData.adaptData{1, i};
        load([scriptDir '/RegressionAnalysis/RegModelResults_',nw ,'/IndvResults/',subID{i} 'defaultsplit_1models_ver00.mat'])
        [~,~,~,Data{i},~] = adaptDataSubject.plotCheckerboards(newLabelPrefix,refEpPost1Early,fh,ph(1,i),refEp,flip); %|EMG_earlyPost1 -  EMG_Baseline
        title({[adaptDataSubject.subData.ID] ['Norm=', num2str(norm(reshape(Data{i},[],1))), '| \beta_{adapt}=', num2str(fitTrans1NoConst.Coefficients.Estimate(1))]});
        
        
    end
    summFlag='nanmedian';
    [~,~,~,Data{i+1},~] = normalizedGroupData.plotCheckerboards(newLabelPrefix,refEpPost1Early,fh,ph(1,n_subjects+1),refEp,flip,summFlag); %|EMG_earlyPost1 -  EMG_Baseline   
    Data{i+1} = nanmedian(Data{i+1}, 4);
    load([scriptDir '/RegressionAnalysis/RegModelResults_',nw ,'/GroupResults/', groupID,'defaultsplit_1flip_1_group_models_ver00.mat'])
    title({['Group'] ['Norm=', num2str(norm(reshape(Data{i+1},[],1))),'| \beta_{adapt}=', num2str(fitTrans1NoConst.Coefficients.Estimate(1))]});
end

set(ph(:,1),'CLim',[-1 1]*1);
set(ph(:,2:end),'YTickLabels',{},'CLim',[-1 1]*1);
set(ph,'FontSize',8)
pos=get(ph(1,end),'Position');
axes(ph(1,end))
colorbar
set(ph(1,end),'Position',pos);
set(gcf,'color','w');


%% Comparing baseline late to Post 1 and Post 2 late

OGbase=defineReferenceEpoch('OG base',ep);
TMbase= defineReferenceEpoch('TM base',ep);
Pos1_Late=defineReferenceEpoch('Post1_{Late}',ep);
Pos2_Late=defineReferenceEpoch('Post2_{Late}',ep);

EpochsOfInteres={OGbase,TMbase,Pos1_Late,Pos2_Late};


if plotIndSubjects
    plotEpochsPlusNorm(EpochsOfInteres,normalizedGroupData,newLabelPrefix,1,0)
end

if plotGroup
    plotEpochsPlusNorm(EpochsOfInteres,normalizedGroupData,newLabelPrefix,0,1)
    
end





