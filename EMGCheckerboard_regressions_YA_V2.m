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
groupID = 'ATR';
saveResAndFigure = true;
plotAllEpoch = true;
plotIndSubjects = true;
plotGroup = true;


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
refEpSlow = defineReferenceEpoch('TM slow',ep);
refEp2 = defineReferenceEpoch('OG base',ep);

newLabelPrefix = defineMuscleList(muscleOrder);
normalizedGroupData = GroupData.normalizeToBaselineEpoch(newLabelPrefix,refEp2); %Normalized by TM base (aka mid baseline)

ll=normalizedGroupData.adaptData{1}.data.getLabelsThatMatch('^Norm');


l2=regexprep(regexprep(ll,'^Norm',''),'_s','');
normalizedGroupData=normalizedGroupData.renameParams(ll,l2);

newLabelPrefix = regexprep(newLabelPrefix,'_s','');

%% Removal of muscle with porblems during data collection
% Remove RF and HIP
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
if plotGroup
    %     ep=defineEpochs_regressionYA('nanmedian');
    summaryflag='nanmedian';
else
    summaryflag=[];
end
if plotAllEpoch
    for i = 1:n_subjects
        
        if plotGroup
            adaptDataSubject = normalizedGroupData;
            figureSaveId = groupID;
        else
            adaptDataSubject = normalizedGroupData.adaptData{1, i};
            figureSaveId = subID{i};
        end
        
        %         adaptDataSubject = normalizedGroupData.adaptData{1, i};
        
        fh=figure('Units','Normalized','OuterPosition',[0 0 1 1]);
        ph=tight_subplot(1,length(ep)-1,[.03 .005],.04,.04);
        flip=true;
        
        
        %         adaptDataSubject.plotCheckerboards(newLabelPrefix,refEp,fh,ph(1,1),[],flip); %plot TM base reference
        %         adaptDataSubject.plotCheckerboards(newLabelPrefix,ep(3,:),fh,ph(1,2),[],flip); %plot TR base reference
        refEp = [];
        adaptDataSubject.plotCheckerboards(newLabelPrefix,ep(1:8,:),fh,ph(1,1:8),refEp,flip,summaryflag);%plot all epochs normalized by the mid baseline
        adaptDataSubject.plotCheckerboards(newLabelPrefix,ep(10:17,:),fh,ph(1,9:length(ep)-1),refEp,flip,summaryflag);%plot all epochs normalized by the mid baseline
        
        %         adaptDataSubject.plotCheckerboards(newLabelPrefix,ep(9,:),fh,ph(1,10),refEpLate,flip);%plot the early Post - late Ada block
        %         title('Post1_{late}-Adapt_{SS}')
        %         adaptDataSubject.plotCheckerboards(newLabelPrefix,ep(10,:),fh,ph(1,11),refEp,flip);%plot all remaining epochs normalized by the fast baseline
        %         title('ShortPos_{early}')
        %
        %         adaptDataSubject.plotCheckerboards(newLabelPrefix,ep(11,:),fh,ph(1,12),refEp,flip);%plot all remaining epochs normalized by the fast baseline
        %         title('NegPos_{early}')
        
        
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

% clc;
% close all
usefft = 0; %Turn on to use the transpose of the EMG(+) as the estimated for EMG(-)
normalizeData = 0; %Vector lenght normalization

%type of regression that we want to run
regModelVersion='default';
% regModelVersion='Adaptive_Feedback'
% regModelVersion='Adaptive_EnvTransition';

% In case that we have multiple cases. Let see if this needs to be updated
splitCount=1;

%For group data try to use the median to be consistnace with Pablo's

summaryflag='nanmean';
%  summaryflag='nanmedian';

% Defining the epochs of interest
ep=defineEpochs_regressionYA(summaryflag);

%Getting the names for specific epocjhs
OGbase=defineReferenceEpoch('OG base',ep);
TMbase= defineReferenceEpoch('TM base',ep); %fast tied 1 if short split 1, slow tied if 2nd split
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
    %         close all;
    for i = 1:n_subjects
        
        
        
        [Data,regressorNames,fh]=RegressionsGeneralization(newLabelPrefix,normalizedGroupData,[],1,0,NegShort,TMbeforeNeg,PosShort,TMbeforePos,...
            AdaptLate,Pos1_Early,Pos1_Late,Pos2_Early, AfterPos, OGbase, TMbase, i);
        
        nw=datestr(now,'yy-mm-dd');
        
        resDir = [scriptDir '/RegressionAnalysis/RegModelResults_', nw,'/IndvResults/'];
        if (saveResAndFigure)
            if not(isfolder(resDir))
                mkdir(resDir)
            end
            saveas(fh, [resDir subID{i} '_Checkerboard_ver' num2str(usefft) num2str(normalizeData) regModelVersion 'split_' num2str(splitCount) '.png'])
            saveas(fh, [resDir subID{i} '_Checkerboard_ver' num2str(usefft) num2str(normalizeData) regModelVersion 'split_' num2str(splitCount)],'epsc')
        end
        
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
%% plot checkerboards and run regression per group
if length(subID) > 1 || plotGroup
    
    
    flip = 1;
      [Data,regressorNames]=RegressionsGeneralization(newLabelPrefix,normalizedGroupData,[],0,1,NegShort,TMbeforeNeg,PosShort,TMbeforePos,...
        AdaptLate,Pos1_Early,Pos1_Late,Pos2_Early, AfterPos, OGbase, TMbase,[]);

        
        nw=datestr(now,'yy-mm-dd');
        resDir = [scriptDir '/RegressionAnalysis/RegModelResults_',nw ,'/GroupResults/'];
        if (saveResAndFigure)
            if not(isfolder(resDir))
                mkdir(resDir)
            end
            saveas(fh, [resDir groupID '_group_Checkerboard_ver' num2str(usefft) num2str(normalizeData) regModelVersion 'split_' num2str(splitCount) 'flip_' num2str(flip) '.png'])
            saveas(fh, [resDir groupID '_group_Checkerboard_ver' num2str(usefft) num2str(normalizeData) regModelVersion 'split_' num2str(splitCount) 'flip_' num2str(flip)], 'epsc')
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

if plotIndSubjects
    for i = 1:n_subjects
        adaptDataSubject = normalizedGroupData.adaptData{1, i};
        load([scriptDir '/RegressionAnalysis/RegModelResults_',nw ,'/IndvResults/',subID{i} 'defaultsplit_1models_ver00.mat'])
        [~,~,~,Data{i},~] = adaptDataSubject.plotCheckerboards(newLabelPrefix,refEpPost1Early,fh,ph(1,i),refEp,flip); %|EMG_earlyPost1 -  EMG_Baseline
        if contains(groupID,'NTR')
            title({[adaptDataSubject.subData.ID] ['Norm=', num2str(norm(reshape(Data{i},[],1))), '| \beta_{adapt}=', num2str(fitTrans1NoConst.Coefficients.Estimate(1))]});
        else
            title({[adaptDataSubject.subData.ID] ['Norm=', num2str(norm(reshape(Data{i},[],1))),'| \beta_{adapt}=', num2str(fitTrans1NoConst.Coefficients.Estimate(1))]});
        end
        
    end
    summFlag='nanmedian';
    [~,~,~,Data{i+1},~] = normalizedGroupData.plotCheckerboards(newLabelPrefix,refEpPost1Early,fh,ph(1,n_subjects+1),refEp,flip,summFlag); %|EMG_earlyPost1 -  EMG_Baseline
    %     vec_norm = norm(Data{i+1});
    %     summFlag='nanmedian';
    %     eval(['fun=@(x) ' summFlag '(x,4);']);
    %     Data{i+1}=fun(Data{i+1});
    Data{i+1} = nanmedian(Data{i+1}, 4);
    load([scriptDir '/RegressionAnalysis/RegModelResults_',nw ,'/GroupResults/', groupID,'defaultsplit_1flip_1_group_models_ver00.mat'])
    if contains(groupID,'NTR')
        title({['Group'] ['Norm=', num2str(norm(reshape(Data{i+1},[],1))),'| \beta_{adapt}=', num2str(fitTrans1NoConst.Coefficients.Estimate(1))]});
    else
        title({['Group'] ['Norm=', num2str(norm(reshape(Data{i+1},[],1))),'| \beta_{adapt}=', num2str(fitTrans1NoConst.Coefficients.Estimate(1))]});
    end
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

if plotIndSubjects
    for i = 1:n_subjects
        adaptDataSubject = normalizedGroupData.adaptData{1, i};
        
        fh=figure('Units','Normalized','OuterPosition',[0 0 1 1],'NumberTitle', 'off', 'Name',subID{i});
        ph=tight_subplot(1,4,[.03 .005],.04,.04);
        flip=true;
        
        [~,~,~,Data{1},~] = adaptDataSubject.plotCheckerboards(newLabelPrefix,OGbase,fh,ph(1,1),[],flip);
        vec_norm = norm(Data{1});
        title({[OGbase.Condition{1}] ['Norm=', num2str(norm(reshape(Data{1},[],1)))]});
        [~,~,~,Data{2},~] = adaptDataSubject.plotCheckerboards(newLabelPrefix,TMbase,fh,ph(1,2),[],flip);
        vec_norm = norm(Data{2});
        title({[TMbase.Condition{1}] ['Norm=', num2str(norm(reshape(Data{2},[],1)))]});
        [~,~,~,Data{3},~] = adaptDataSubject.plotCheckerboards(newLabelPrefix,Pos1_Late,fh,ph(1,3),[],flip);
        vec_norm = norm(Data{3});
        title({[Pos1_Late.Condition{1}] ['Norm=', num2str(norm(reshape(Data{3},[],1)))]});
        [~,~,~,Data{4},~] = adaptDataSubject.plotCheckerboards(newLabelPrefix,Pos2_Late,fh,ph(1,4),[],flip);
        vec_norm = norm(Data{4});
        title({[Pos2_Late.Condition{1}] ['Norm=', num2str(norm(reshape(Data{4},[],1)))]});
        
        set(ph(:,1),'CLim',[-1 1]*1);
        set(ph(:,2:end),'YTickLabels',{},'CLim',[-1 1]*1);
        set(ph,'FontSize',8)
        pos=get(ph(1,end),'Position');
        axes(ph(1,end))
        colorbar
        set(ph(1,end),'Position',pos);
        set(gcf,'color','w');
        
    end
    
    
end

if plotGroup
    
    fh=figure('Units','Normalized','OuterPosition',[0 0 1 1]);
    ph=tight_subplot(1,4,[.03 .005],.04,.04);
    flip=true;
    summFlag='nanmedian';
    %     eval(['fun=@(x) ' summFlag '(x,4);']);
    
    [~,~,~,Data{1},~] = normalizedGroupData.plotCheckerboards(newLabelPrefix,OGbase,fh,ph(1,1),[],flip,summFlag);
    %     Data{1}=fun(Data{1});
    Data{1} = nanmedian(Data{1}, 4);
    title({[OGbase.Condition{1}] ['Norm=', num2str(norm(reshape(Data{1},[],1)))]});
    
    [~,~,~,Data{2},~] = normalizedGroupData.plotCheckerboards(newLabelPrefix,TMbase,fh,ph(1,2),[],flip,summFlag);
    %     Data{2}=fun(Data{2});
    Data{2} = nanmedian(Data{2}, 4);
    title({[TMbase.Condition{1}] ['Norm=', num2str(norm(reshape(Data{2},[],1)))]});
    
    [~,~,~,Data{3},~] = normalizedGroupData.plotCheckerboards(newLabelPrefix,Pos1_Late,fh,ph(1,3),[],flip,summFlag);
    %     Data{3}=fun(Data{3});
    Data{3} = nanmedian(Data{3}, 4);
    title({[Pos1_Late.Condition{1}] ['Norm=', num2str(norm(reshape(Data{3},[],1)))]});
    
    [~,~,~,Data{4},~] = normalizedGroupData.plotCheckerboards(newLabelPrefix,Pos2_Late,fh,ph(1,4),[],flip,summFlag);
    %     Data{4}=fun(Data{4});
    Data{4} = nanmedian(Data{4}, 4);
    title({[Pos2_Late.Condition{1}] ['Norm=', num2str(norm(reshape(Data{4},[],1)))]});
    
    
    
end
set(ph(:,1),'CLim',[-1 1]*1);
set(ph(:,2:end),'YTickLabels',{},'CLim',[-1 1]*1);
set(ph,'FontSize',8)
pos=get(ph(1,end),'Position');
axes(ph(1,end))
colorbar
set(ph(1,end),'Position',pos);
set(gcf,'color','w');





%%
%    vec_norm = vecnorm(Data{1});
%% plot subsets of epochs: AE with context specific baseline correction
%AE only pertains to session 1 and long protocols.

% refEpOG = defineReferenceEpoch('OGbase',ep);
% refEpTR = defineReferenceEpoch('TRbase',ep);
% post1ep = ep(strcmp(ep.Properties.ObsNames,'Post1_{Early}'),:);
% post2ep = ep(strcmp(ep.Properties.ObsNames,'Post2_{Early}'),:);
% for flip = [1,2] %2 legs separate first (flip = 1) and then asymmetry (flip = 2)
%     %     the flip asymmetry plots average of summation and the average of
%     %     asymmetry.
%     for i = 1:n_subjects
%         if plotGroup
%             adaptDataSubject = normalizedTMFullAbrupt;
%             figureSaveId = groupID;
%         else
%             adaptDataSubject = normalizedTMFullAbrupt.adaptData{1, i};
%             figureSaveId = subID{i};
%         end
%
%         fh=figure('Units','Normalized','OuterPosition',[0 0 1 1]);
%         ph=tight_subplot(1,4,[.03 .005],.04,.04);
%
%         % plot after effects only
%         adaptDataSubject.plotCheckerboards(newLabelPrefix,refEpOG,fh,ph(1,1),[],flip); %plot OG base
%         adaptDataSubject.plotCheckerboards(newLabelPrefix,refEpTR,fh,ph(1,2),[],flip); %plot OG base with shoe
%
%         if (contains(groupID, 'TS'))%correct post 1 with OG no nimbus, post 2 with OG with nimbus, i.e.,TR
%             adaptDataSubject.plotCheckerboards(newLabelPrefix,post1ep,fh,ph(1,3),refEpOG,flip); %post1 is OG
%             adaptDataSubject.plotCheckerboards(newLabelPrefix,post2ep,fh,ph(1,4),refEpTR,flip); %post2 is with Nimbus(TR)
%         elseif (contains(groupID, 'TR')) %correct post 1 with TR(nimbus), post 2 with OG (No nimbus)
%             adaptDataSubject.plotCheckerboards(newLabelPrefix,post1ep,fh,ph(1,3),refEpTR,flip);
%             adaptDataSubject.plotCheckerboards(newLabelPrefix,post2ep,fh,ph(1,4),refEpOG,flip);
%         end
%
%         set(ph(:,1),'CLim',[-1 1]);
%         set(ph(:,2:end),'YTickLabels',{},'CLim',[-1 1]*0.5);
%         set(ph,'FontSize',8)
%         pos=get(ph(1,end),'Position');
%         axes(ph(1,end))
%         colorbar
%         set(ph(1,end),'Position',pos);
%         set(gcf,'color','w');
%
%         if (saveResAndFigure)
%             resDir = [scriptDir '/RegressionAnalysis/RegModelResults_V14/'];
%             if plotGroup
%                 resDir = [resDir 'GroupResults/'];
%             end
%             if not(isfolder(resDir))
%                 mkdir(resDir)
%             end
%             if flip == 1
%                 saveas(fh, [resDir figureSaveId '_CheckerBoard_AE_SpecificBase.png'])
%                 %                 saveas(fh, [resDir figureSaveId '_AllEpochCheckerBoard_' num2str(session)],'epsc')
%             else
%                 saveas(fh, [resDir figureSaveId '_CheckerBoard_AE_SpecificBase_Asym.png'])
%             end
%         end
%
%         if plotGroup
%             break
%         end
%     end
% end

%% plot subsets of epochs: AE with context specific baseline correction
%AE only pertains to session 1 and long protocols.
% refEpOG = defineReferenceEpoch('OGbase',ep);
% refEpTR = defineReferenceEpoch('TM base',ep);
% post1ep = ep(strcmp(ep.Properties.ObsNames,'Post1_{Early}'),:);
% post2ep = ep(strcmp(ep.Properties.ObsNames,'Post2_{Early}'),:);
% post1lateep = ep(strcmp(ep.Properties.ObsNames,'Post1_{Late}'),:);
% post2lateep = ep(strcmp(ep.Properties.ObsNames,'Post2_{Late}'),:);
%
% Data = cell(1,4);
%
% for flip = [1]%,2] %2 legs separate first (flip = 1) and then asymmetry (flip = 2)
%     %     the flip asymmetry plots average of summation and the average of
%     %     asymmetry.
%     for i = 1:n_subjects
%         if plotGroup
%             adaptDataSubject = normalizedTMFullAbrupt;
%             figureSaveId = groupID;
%         else
%             adaptDataSubject = normalizedTMFullAbrupt.adaptData{1, i};
%             figureSaveId = subID{i};
%         end
%
%         fh=figure('Units','Normalized','OuterPosition',[0 0 1 1]);
%         ph=tight_subplot(1,8)%,[.03 .005],.04,.04);
%
%         % plot after effects only
%         [~,~,~,Data{1},~]= adaptDataSubject.plotCheckerboards(newLabelPrefix,refEpOG,fh,ph(1,1),[],flip); %plot OG base
%         [~,~,~,Data{2},~]= adaptDataSubject.plotCheckerboards(newLabelPrefix,refEpTR,fh,ph(1,2),[],flip); %plot OG base with shoe
%
%         %         if (contains(groupID, 'TS'))%correct post 1 with OG no nimbus, post 2 with OG with nimbus, i.e.,TR
%         %             adaptDataSubject.plotCheckerboards(newLabelPrefix,post1ep,fh,ph(1,3),refEpOG,flip); %post1 is OG
%         %             adaptDataSubject.plotCheckerboards(newLabelPrefix,post2ep,fh,ph(1,4),refEpTR,flip); %post2 is with Nimbus(TR)
%         %         elseif (contains(groupID, 'TR')) %correct post 1 with TR(nimbus), post 2 with OG (No nimbus)
%         %             adaptDataSubject.plotCheckerboards(newLabelPrefix,post1ep,fh,ph(1,3),refEpTR,flip);
%         %             adaptDataSubject.plotCheckerboards(newLabelPrefix,post2ep,fh,ph(1,4),refEpOG,flip);
%         %         end
%
%         [~,~,~,Data{3},~]= adaptDataSubject.plotCheckerboards(newLabelPrefix,post1ep,fh,ph(1,3),[],flip);
%         [~,~,~,Data{4},~]= adaptDataSubject.plotCheckerboards(newLabelPrefix,post2ep,fh,ph(1,4),[],flip);
%         [~,~,~,Data{5},~]=adaptDataSubject.plotCheckerboards(newLabelPrefix,post1lateep,fh,ph(1,5),[],flip);
%         [~,~,~,Data{6},~]=adaptDataSubject.plotCheckerboards(newLabelPrefix,post2lateep,fh,ph(1,6),[],flip);
%
%         %         if strcmpi(regModelVersion, 'TR')
%         adaptDataSubject.plotCheckerboards(newLabelPrefix,post1lateep,fh,ph(1,7),refEpTR,flip);
%         title('Post1_{Late}-TRbase_{late}')
%
%         disp('Cosine(Post1_{Late},TRbase_{late})')
%         cos= Cosine2Matrix(Data{5},Data{2})
%
%
%
%         adaptDataSubject.plotCheckerboards(newLabelPrefix,post2lateep,fh,ph(1,8),refEpOG,flip);
%         title('Post2_{Late}-OGbase_{late}')
%
%
%         disp('Cosine(Post2_{Late},OGbase_{late})')
%         cos= Cosine2Matrix(Data{1},Data{6})
%
%         %         else
%         %
%         %             adaptDataSubject.plotCheckerboards(newLabelPrefix,post1lateep,fh,ph(1,7),refEpOG,flip);
%         %             title('Post1_{Late}-OGbase_{late}')
%         %             disp('Cosine(Post1_{Late},OGbase_{late})')
%         %             cos= Cosine2Matrix(Data{1},Data{5})
%         %
%         %             %             cosine(Data{1},Data{3})
%         %
%         %             adaptDataSubject.plotCheckerboards(newLabelPrefix,post2lateep,fh,ph(1,8),refEpTR,flip);
%         %             title('Post2_{Late}-TRbase_{late}')
%         %             disp('Cosine(Post2_{Late},TRbase_{late})')
%         %             cos= Cosine2Matrix(Data{2},Data{6})
%         %             %             cosine(Data{2},Data{4})
%         %         end
%
%
%         %         adaptDataSubject.plotCheckerboards(newLabelPrefix,post1ep,fh,ph(1,7),ep(4,:),flip);
%         %         title('trans2');
%         %         adaptDataSubject.plotCheckerboards(newLabelPrefix,post2ep,fh,ph(1,8),post1lateep,flip);
%         %         title('trans2');
%         %         adaptDataSubject.plotCheckerboards(newLabelPrefix,refEpOG,fh,ph(1,9),refEpTR,flip);
%         %         title('OG-TRbase');
%
%         set(ph(:,1),'CLim',[-1 1]);
%         set(ph(:,2:end),'YTickLabels',{},'CLim',[-1 1]);
%         set(ph,'FontSize',8)
%         pos=get(ph(1,end),'Position');
%         axes(ph(1,end))
%         colorbar
%         set(ph(1,end),'Position',pos);
%         set(gcf,'color','w');
%
%         currLabelPrefix = newLabelPrefix;
%
%
%
%         if (saveResAndFigure)
%             resDir = [scriptDir '/RegressionAnalysis/RegModelResults_V15/'];
%             if plotGroup
%                 resDir = [resDir 'GroupResults/'];
%             end
%             if not(isfolder(resDir))
%                 mkdir(resDir)
%             end
%             if flip == 1
%                 saveas(fh, [resDir figureSaveId '_CheckerBoard_AE_NoBase.png'])
%                 %                 saveas(fh, [resDir figureSaveId '_AllEpochCheckerBoard_' num2str(session)],'epsc')
%             else
%                 saveas(fh, [resDir figureSaveId '_CheckerBoard_AE_NoBase_Asym.png'])
%             end
%         end
%
%         if plotGroup
%             break
%         end
%     end
% end
%
%
%
%% Compare Off perturbation to OG vs  off pertubation + multiEnvSwitch
% Compare TRSplitLate - OGPostEarly ~ (TRBase - OGBase) + (TRSplitLate - TRPostEarly)
% ep = defineEpochVR_OG_UpdateV3('nanmean', subID);
% %         close all;
% for i = 1:n_subjects
%     adaptDataSubject = normalizedTMFullAbrupt.adaptData{1, i};
%
%     fh=figure('Units','Normalized','OuterPosition',[0 0 1 1]);
%     ph=tight_subplot(1,3,[.03 .005],.04,.04);
%     flip=true;
%
%     Data = {};
%     %all labels should be the same, no need to save again.
%     %  OGPost - TRSplitLate
%     [~,~,labels,Data{1},dataRef2]=adaptDataSubject.plotCheckerboards(newLabelPrefix,ep(15,:),fh,ph(1,1),ep(14,:),flip);
%     title('OGPost - PosSplitLate(fastPrior)')
%     % OGBase - TRBase
%     [~,~,~,Data{2},~] = adaptDataSubject.plotCheckerboards(newLabelPrefix,ep(1,:),fh,ph(1,2),ep(2,:),flip); %  -(TR base - OG base) = OG base - TR base, env switching
%     title('OGBaseLate - TRBaseLate')
%     % TRBase - TRSplitEarly
%     [~,~,~,Data{3},~]=adaptDataSubject.plotCheckerboards(newLabelPrefix,ep(2,:),fh,ph(1,3),ep(4,:),flip);
%     title('TMBaseLate - TMSplitEarly')
%
%     set(ph(:,1),'CLim',[-1 1]);
%     set(ph(:,2:end),'YTickLabels',{},'CLim',[-1 1]*1.5);
%     set(ph,'FontSize',8)
%     pos=get(ph(1,end),'Position');
%     axes(ph(1,end))
%     colorbar
%     set(ph(1,end),'Position',pos);
%     set(gcf,'color','w');
%
%     resDir = [scriptDir '/RegressionAnalysis/RegModelResults_V14/'];
%     if saveResAndFigure
%         if not(isfolder(resDir))
%             mkdir(resDir)
%         end
%         saveas(fh, [resDir subID{i} '_TROffLinear' '.png'])
% %             saveas(fh, [resDir subID{i} '_ShoeOffLinear'],'epsc')
%     end
%
%     % run regression and save results
%     format compact % format loose %(default)
%     YvsTerm1Correlation = corrcoef(Data{1},Data{2})
%     YvsTerm2Correlation = corrcoef(Data{1},Data{3})
%     corr_coef={YvsTerm1Correlation, YvsTerm2Correlation};
%
%     for j = 1:size(Data,2)
%         Data{j} = reshape(Data{j}, [],1); %make it a column vector
%     end
%     YvsTerm1Cos = cosine(Data{1},Data{2})
%     YvsTerm2Cos = cosine(Data{1},Data{3})
%     cosine_values={YvsTerm1Cos, YvsTerm2Cos};
%
%     %%% Run regression to see if the LHS and RHS are equal
%     tableData=table(Data{1},Data{2},Data{3},'VariableNames',{'OGToTRSplit', 'OGToTR', 'TRToSplit'});
%     linearityAssessmentModel=fitlm(tableData, 'OGToTRSplit~OGToTR+TRToSplit-1')%exclude constant
%     Rsquared = linearityAssessmentModel.Rsquared
%
%     if saveResAndFigure
%         save([resDir subID{i} '_TROffLinear_CorrCoef_Model'],'corr_coef','cosine_values','linearityAssessmentModel')
%     end
% end