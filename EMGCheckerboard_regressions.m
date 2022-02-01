%This file will run EMG regressions and plot checkerboards.
% Version update from V3: focus on subjects that performed 2 sets of
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
groupID = 'C';
saveResAndFigure = false;
plotAllEpoch = true;
plotIndSubjects = true;
plotGroup = true;
kinenatics=true;

scriptDir = fileparts(matlab.desktop.editor.getActiveFilename);
if kinenatics==1
    files = dir ([scriptDir '/data_Kinematics/' groupID '*params.mat']);
    
else
    
    files = dir ([scriptDir '/data_reprocess2021/' groupID '*params2021.mat']);
end
n_subjects = size(files,1);
subID = cell(1, n_subjects);
sub=cell(1,n_subjects);
for i = 1:n_subjects
    sub{i} = files(i).name;
    if kinenatics==1
        subID{i} = sub{i}(1:end-10);
    else
        subID{i} = sub{i}(1:end-14);
    end
end
subID

regModelVersion =  'default'; %'default'
% if (contains(groupID, 'C') || contains(groupID, 'VROG'))
%     regModelVersion = 'TS'
% elseif (contains(groupID, 'CTR'))
% regModelVersion = 'TR'
% end

%% load and prep data
normalizedTMFullAbrupt=adaptationData.createGroupAdaptData(sub);
normalizedTMFullAbrupt=normalizedTMFullAbrupt.removeBadStrides;

muscleOrder={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT', 'VM', 'VL', 'RF', 'TFL', 'GLU', 'HIP','ADM'};

n_muscles = length(muscleOrder);

n_subjects = length(subID);

ep=defineEpochs_regression('nanmean');
refEp = defineReferenceEpoch('TM base',ep);
refEpLate = defineReferenceEpoch('Adaptation',ep);
refEpSlow = defineReferenceEpoch('TM slow',ep);
refEp2 = defineReferenceEpoch('OGbase',ep);

newLabelPrefix = defineMuscleList(muscleOrder);


normalizedTMFullAbrupt = normalizedTMFullAbrupt.normalizeToBaselineEpoch(newLabelPrefix,refEp); %Normalized by TM base (aka mid baseline)

ll=normalizedTMFullAbrupt.adaptData{1}.data.getLabelsThatMatch('^Norm');
l2=regexprep(regexprep(ll,'^Norm',''),'_s','s');
normalizedTMFullAbrupt=normalizedTMFullAbrupt.renameParams(ll,l2);
newLabelPrefix = regexprep(newLabelPrefix,'_s','s');

%% EMG aftereffects

refEpPost1Early= defineReferenceEpoch('Post1_{Early}',ep);
if contains(groupID,'C')
    refEp= defineReferenceEpoch('TM base',ep); %fast tied 1 if short split 1, slow tied if 2nd split
end

fh=figure('Units','Normalized','OuterPosition',[0 0 1 1]);
ph=tight_subplot(1,n_subjects+1,[.03 .005],.04,.04);
flip = [2];

% nw=datestr(now,'yy-mm-dd');
% nw='22-01-31';
if plotIndSubjects
    for i = 1:n_subjects
        adaptDataSubject = normalizedTMFullAbrupt.adaptData{1, i};
        [~,~,~,Data{i},~] = adaptDataSubject.plotCheckerboards(newLabelPrefix,refEpPost1Early,fh,ph(1,i),refEp,flip); %|EMG_earlyPost1 -  EMG_Baseline
        title({[adaptDataSubject.subData.ID] ['Norm=', num2str(norm(reshape(Data{i},[],1)))]});
        
        
    end
    summFlag='nanmedian';
    [~,~,~,Data{i+1},~] = normalizedTMFullAbrupt.plotCheckerboards(newLabelPrefix,refEpPost1Early,fh,ph(1,n_subjects+1),refEp,flip,summFlag); %|EMG_earlyPost1 -  EMG_Baseline   
    Data{i+1} = nanmedian(Data{i+1}, 4);
%     load([scriptDir '/RegressionAnalysis/RegModelResults_',nw ,'/GroupResults/', groupID,'defaultsplit_1flip_1_group_models_ver00.mat'])
    title({['Group'] ['Norm=', num2str(norm(reshape(Data{i+1},[],1)))]});
end

set(ph(:,1),'CLim',[-1 1]*1);
set(ph(:,2:end),'YTickLabels',{},'CLim',[-1 1]*1);
set(ph,'FontSize',8)
pos=get(ph(1,end),'Position');
axes(ph(1,end))
colorbar
set(ph(1,end),'Position',pos);
set(gcf,'color','w');

demographics= normalizedTMFullAbrupt.GroupDemographics;




%%
% Remove fRF from CTS
% if contains(groupID,'CTS')
% %     badMuscleNames = {'fRFs','sTAs'};
%     badMuscleNames = {'sTAs'};
%     badMuscleIdx=[];
%     for bm = badMuscleNames
%         badMuscleIdx = [badMuscleIdx, find(ismember(newLabelPrefix,bm))];
%     end
%     newLabelPrefix = newLabelPrefix(setdiff(1:end, badMuscleIdx))
% elseif contains(groupID,'CTR_02')
%     badMuscleNames = {'fTFLs'};
%     badMuscleIdx=[];
%     for bm = badMuscleNames
%         badMuscleIdx = [badMuscleIdx, find(ismember(newLabelPrefix,bm))];
%     end
%     newLabelPrefix = newLabelPrefix(setdiff(1:end, badMuscleIdx))
% elseif contains(groupID,'CTR_05')
%         badMuscleNames = {'sHIPs'};
%         badMuscleIdx=[];
%         for bm = badMuscleNames
%             badMuscleIdx = [badMuscleIdx, find(ismember(newLabelPrefix,bm))];
%         end
%     newLabelPrefix = newLabelPrefix(setdiff(1:end, badMuscleIdx))
% elseif contains(groupID,'CTR')
%         badMuscleNames = {'sHIPs','fTFLs','fHIPs','sTFLs'};
% %          badMuscleNames = {'fTFLs','sRFs'};
%         badMuscleIdx=[];
%         for bm = badMuscleNames
%             badMuscleIdx = [badMuscleIdx, find(ismember(newLabelPrefix,bm))];
%         end
%         newLabelPrefix = newLabelPrefix(setdiff(1:end, badMuscleIdx))
% end
%%
% adaptDataSubject = normalizedTMFullAbrupt.adaptData{1, 1};
%
% fh=figure('Units','Normalized','OuterPosition',[0 0 1 1]);
% ph=tight_subplot(1,1,[.03 .005],.04,.04);
% flip=true;
%
% adaptDataSubject.plotCheckerboards(newLabelPrefix,defineReferenceEpoch('OGbase',ep),fh,ph(1,1),[],flip); %plot TM slow
% % adaptDataSubject.plotCheckerboards(newLabelPrefix,refEp,fh,ph(1,2),[],flip); %plot TR base reference
%
%
% % adaptDataSubject.plotCheckerboards(newLabelPrefix,defineReferenceEpoch('TMslow',ep),fh,ph(1,1),[],flip); %plot TR base reference
%
%% plot epochs
normalize=0;
if plotAllEpoch
    for i = 1%:n_subjects
        
        if plotGroup
            adaptDataSubject = normalizedTMFullAbrupt;
            figureSaveId = groupID;
        else
            adaptDataSubject = normalizedTMFullAbrupt.adaptData{1, i};
            figureSaveId = subID{i};
        end
        
        adaptDataSubject = normalizedTMFullAbrupt.adaptData{1, i};
        
        fh=figure('Units','Normalized','OuterPosition',[0 0 1 1]);
        ph=tight_subplot(1,length(ep)+1,[.03 .005],.04,.04);
        flip=true;
        
        
        adaptDataSubject.plotCheckerboards(newLabelPrefix,refEp,fh,ph(1,1),[],flip); %plot TM base reference
        %         adaptDataSubject.plotCheckerboards(newLabelPrefix,ep(3,:),fh,ph(1,2),[],flip); %plot TR base reference
        
        adaptDataSubject.plotCheckerboards(newLabelPrefix,ep(1:8,:),fh,ph(1,2:9),refEp,flip);%plot all epochs normalized by the mid baseline
        adaptDataSubject.plotCheckerboards(newLabelPrefix,ep(9,:),fh,ph(1,10),refEpLate,flip);%plot the early Post - late Ada block
        title('Post1_{late}-Adapt_{SS}')
        adaptDataSubject.plotCheckerboards(newLabelPrefix,ep(10,:),fh,ph(1,11),refEp,flip);%plot all remaining epochs normalized by the fast baseline
        title('ShortPos_{early}-TMbase_{late}')
        
        
        set(ph(:,1),'CLim',[-1 1]*1);
        %     set(ph(:,2:end),'YTickLabels',{},'CLim',[-1 1]*2);
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

%% Regressor V2, prepare data for regressor checkerboards and regression model
% - Update from disucssion on Tuesday MAy 4, 2021
% - Refer to OneNote for naming detailsxc

% Adapt VR- baseline VR
% OG base - baseline VR
% baseline - EMG_split(-)
clc;
close all
usefft = 1; normalizeData = 0;
normalize=0;
regModelVersion='default';
% regModelVersion='Adaptive_EnvTransition';
splitCount=1
% for splitCount = 1
%     splitCount
%     if splitCount == 1 %first short split, fast baseline - split
%         ep=defineEpochVR_OG_UpdateV3('nanmean');
%         refPosShort = defineReferenceEpoch('WithinEnvSwitch (-\DeltaEMG_{on(+)})', ep);
%     elseif splitCount == 2%2nd short split, slow baseline - split
%         ep=defineEpochVR_OG_UpdateV4('nanmean');
%         refPosShort = defineReferenceEpoch('WithinEnvSwitch (-\DeltaEMG_{on(+)})2', ep);
%     elseif splitCount == 3 %EMG_pos (Pos Short - fast base );  EMG_neg (Neg Short 2 - Slow base )
%         ep=defineEpochVR_OG_UpdateV5('nanmean', groupID);
%         refPosShort = defineReferenceEpoch('WithinEnvSwitch (-\DeltaEMG_{on(+)})', ep);
%     elseif splitCount == 4 %2nd short split positive, slow baseline
%         ep=defineEpochVR_OG_UpdateV6('nanmean');
%         refPosShort = defineReferenceEpoch('WithinEnvSwitch (-\DeltaEMG_{on(+)})', ep);
%     elseif splitCount == 5
%         ep=defineEpochVR_OG_UpdateV7('nanmean');
%         refPosShort = defineReferenceEpoch('WithinEnvSwitch (-\DeltaEMG_{on(+)})', ep);
%     end
summaryflag='nanmean';
ep=defineEpochs_regression(summaryflag);
refEpAdaptLate = defineReferenceEpoch('Adaptation',ep);
refEpPost1Late= defineReferenceEpoch('Post1_{Late}',ep);
refEp= defineReferenceEpoch('TM base',ep); %fast tied 1 if short split 1, slow tied if 2nd split

Rsquared1_Adj=[];
Rsquared1_Ord=[];

Coefficients_trans1=[];
Coefficients_trans2=[];

%% plot checkerboard and run regression per subject
if plotIndSubjects
    %         close all;
    for i = 1:n_subjects
        disp(['Subject=', num2str(i)])
        adaptDataSubject = normalizedTMFullAbrupt.adaptData{1, i};
        
        fh=figure('Units','Normalized','OuterPosition',[0 0 1 1]);
        ph=tight_subplot(1,5,[.03 .005],.04,.04);
        flip=true;
        
        Data = {}; %in order: adapt, dataEnvSwitch, dataTaskSwitch, dataTrans1, dataTrans2
        if usefft
            [~,~,labels,Data{1},dataRef2]=adaptDataSubject.plotCheckerboards(newLabelPrefix,ep(10,:),fh,ph(1,1),refEp,flip); %  EMG_split(-) - TM base
            title('EMG(+)-TM_{base}')
        else
            [~,~,labels,Data{1},dataRef2]=adaptDataSubject.plotCheckerboards(newLabelPrefix,ep(10,:),fh,ph(1,1),refEp,flip); %  EMG_split(-) - TM base
        end
        %all labels should be the same, no need to save again.
        [~,~,~,Data{2},~] = adaptDataSubject.plotCheckerboards(newLabelPrefix,refEp,fh,ph(1,2),ep(10,:),flip); %TM base - EMG_on(+)
        title('TM_{base}-EMG(+)')
        [~,~,~,Data{3},~] = adaptDataSubject.plotCheckerboards(newLabelPrefix,ep(1,:),fh,ph(1,3),refEp,flip); %  OG base - TR base, env switching
        title('OG_{base}-TM_{base}')
        [~,~,~,Data{4},~] = adaptDataSubject.plotCheckerboards(newLabelPrefix,ep(5,:),fh,ph(1,4),refEpAdaptLate,flip); %Post1 - Adaptation_{SS}, transition 1
        title('Post1_{e}-Adapt_{SS}')
        [~,~,~,Data{5},~] = adaptDataSubject.plotCheckerboards(newLabelPrefix,ep(7,:),fh,ph(1,5),refEpPost1Late,flip); %Post2 early - Post 1 late, transition 2
        title('Post2_{e}-Post1{l}')
        %     [~,~,labels,dataE{1},dataRef{1}]=adaptDataSubject.plotCheckerboards(newLabelPrefix,ep,fh,ph(1,2:end),refEp,flip);%Second, the rest:
        
        set(ph(:,1),'CLim',[-1 1]);
        %     set(ph(:,2:end),'YTickLabels',{},'CLim',[-1 1]*2);
        set(ph(:,2:end),'YTickLabels',{},'CLim',[-1 1]*1.5);
        set(ph,'FontSize',8)
        pos=get(ph(1,end),'Position');
        axes(ph(1,end))
        colorbar
        set(ph(1,end),'Position',pos);
        set(gcf,'color','w');
        
        resDir = [scriptDir '/RegressionAnalysis/RegModelResults_V14/'];
        if (saveResAndFigure)
            if not(isfolder(resDir))
                mkdir(resDir)
            end
            %                 saveas(fh, [resDir subID{i} '_Checkerboard_ver' num2str(usefft) num2str(normalizeData) regModelVersion 'split_' num2str(splitCount) '.png'])
            %                 saveas(fh, [resDir subID{i} '_Checkerboard_ver' num2str(usefft) num2str(normalizeData) regModelVersion 'split_' num2str(splitCount)],'epsc')
        end
        
        % run regression and save results
        format compact % format loose %(default)
        % not normalized first, then normalized, arugmnets order: (Data, normalizeData, isGroupData, dataId, resDir, saveResAndFigure, version, usefft)
        fprintf('\n')
        splitCount
        [trans1,trans2]=runRegression_Generalization(Data, false, false, [subID{i} regModelVersion 'split_' num2str(splitCount)], resDir, saveResAndFigure, regModelVersion, usefft)
        runRegression_Generalization(Data, true, false, [subID{i} regModelVersion 'split_' num2str(splitCount)], resDir, saveResAndFigure, regModelVersion, usefft)
Rsquared1_Adj(i,1)=trans1.Rsquared.Adjusted;
Rsquared1_Ord(i,1)=trans1.Rsquared.Ordinary;

Rsquared2_Adj(i,1)=trans2.Rsquared.Adjusted;
Rsquared2_Ord(i,1)=trans2.Rsquared.Ordinary;

Coefficients_trans1(:,i)=trans1.Coefficients.Estimate;
Coefficients_trans2(:,i)=trans2.Coefficients.Estimate;
   
    end

    
    
end
%% plot checkerboards and run regression per group
if length(subID) > 1 || plotGroup
    for flip = [1]
        summaryflag='nanmedian';
        %             if splitCount == 1 %first short split, fast baseline - split
        %                 ep=defineEpochVR_OG_UpdateV3('nanmean');
        %                 refPosShort = defineReferenceEpoch('WithinEnvSwitch (-\DeltaEMG_{on(+)})', ep);
        %             elseif splitCount == 2%2nd short split, slow baseline - split
        %                 ep=defineEpochVR_OG_UpdateV4('nanmean');
        %                 refPosShort = defineReferenceEpoch('WithinEnvSwitch (-\DeltaEMG_{on(+)})2', ep);
        %             elseif splitCount == 3 %EMG_pos (Pos Short - fast base );  EMG_neg (Neg Short 2 - Slow base )
        %                 ep=defineEpochVR_OG_UpdateV5('nanmean', groupID);
        %                 refPosShort = defineReferenceEpoch('WithinEnvSwitch (-\DeltaEMG_{on(+)})', ep);
        %             elseif splitCount == 4  %EMG_pos (Pos Short - fast base );  EMG_neg (Pos Short 2 - Slow base )
        %                 ep=defineEpochVR_OG_UpdateV6('nanmean');
        %                 refPosShort = defineReferenceEpoch('WithinEnvSwitch (-\DeltaEMG_{on(+)})', ep);
        %             elseif splitCount == 5  %EMG_pos (Pos Short - fast base );  EMG_neg (Pos Short  - Slow base )
        %                 ep=defineEpochVR_OG_UpdateV7('nanmean');
        %                 refPosShort = defineReferenceEpoch('WithinEnvSwitch (-\DeltaEMG_{on(+)})', ep);
        %             end
        %             ep=defineEpochs_regression('nanmean');
        %             refEpAdaptLate = defineReferenceEpoch('Adaptation',ep);
        %             refEpPost1Late= defineReferenceEpoch('Post1_{Late}',ep);
        %             refEp= defineReferenceEpoch('TM base',ep); %fast tied 1 if short split 1, slow tied if 2nd split
        
        %             refEpAdaptLate = defineReferenceEpoch('Task_{Switch}',ep);
        %             refEpPost1Late= defineReferenceEpoch('Post1_{Late}',ep);
        %             refEp= defineReferenceEpoch('TMbase',ep); %fast tied 1 if short split 1, slow tied if 2nd split
        
        fh=figure('Units','Normalized','OuterPosition',[0 0 1 1]);
        ph=tight_subplot(1,5,[.03 .005],.04,.04);
        
        Data = {}; %in order: adapt, dataEnvSwitch, dataTaskSwitch, dataTrans1, dataTrans2
        if usefft
            [~,~,labels,Data{1},dataRef2]= normalizedTMFullAbrupt.plotCheckerboards(newLabelPrefix,ep(10,:),fh,ph(1,1),refEp,flip,summaryflag); %  EMG_split(-) - TM base
            title('EMG(+)-TM_{base}')
        else
            [~,~,labels,Data{1},dataRef2]= normalizedTMFullAbrupt.plotCheckerboards(newLabelPrefix,ep(10,:),fh,ph(1,1),refEp,flip,summaryflag); %  EMG_split(-) - TM base
        end
        
        [~,~,~,Data{2},~] = normalizedTMFullAbrupt.plotCheckerboards(newLabelPrefix,refEp,fh,ph(1,2),ep(10,:),flip,summaryflag); %TM base - EMG_on(+)
        title('TM_{base}-EMG(+)')
        
        [~,~,~,Data{3},~] = normalizedTMFullAbrupt.plotCheckerboards(newLabelPrefix,ep(1,:),fh,ph(1,3),refEp,flip,summaryflag); %  OG base - TR base, env switching
        title('OG_{base}-TM_{base}')
        
        [~,~,~,Data{4},~] = normalizedTMFullAbrupt.plotCheckerboards(newLabelPrefix,ep(5,:),fh,ph(1,4),refEpAdaptLate,flip,summaryflag); %Post1 - Adaptation_{SS}, transition 1
        title('Post1_{e}-Adapt_{SS}')
        
        [~,~,~,Data{5},~] = normalizedTMFullAbrupt.plotCheckerboards(newLabelPrefix,ep(7,:),fh,ph(1,5),refEpPost1Late,flip,summaryflag); %Post2 early - Post 1 late, transition 2
        title('Post2_{e}-Post1{l}')
        %all labels should be the same, no need to save again.
        %             [~,~,~,Data{2},~] = normalizedTMFullAbrupt.plotCheckerboards(newLabelPrefix,refPosShort,fh,ph(1,2),ep(4,:),flip); % Noadapt (env-driven), TM base - EMG_on(+)
        %             if ~strcmp(groupID, 'CTR')
        %                 [~,~,~,Data{3},~] = normalizedTMFullAbrupt.plotCheckerboards(newLabelPrefix,ep(1,:),fh,ph(1,3),ep(5,:),flip); % OG base - TR base, env switching
        %             else %CTR TR base is fast but post is mid, so should use mid for env transition
        %                 [~,~,~,Data{3},~] = normalizedTMFullAbrupt.plotCheckerboards(newLabelPrefix,ep(1,:),fh,ph(1,3),refEp,flip); % OG base - TR base, env switching
        %             end
        %             [~,~,~,Data{4},~] = normalizedTMFullAbrupt.plotCheckerboards(newLabelPrefix,ep(8,:),fh,ph(1,4),refEpAdaptLate,flip); %OGafter - Adaptation_{SS}, transition 1
        %             [~,~,~,Data{5},~] = normalizedTMFullAbrupt.plotCheckerboards(newLabelPrefix,ep(11,:),fh,ph(1,5),refEpPost1Late,flip); %TM post VR early - OG post late, transition 2
        %     [~,~,labels,dataE{1},dataRef{1}]=normalizedTMFullAbrupt.plotCheckerboards(newLabelPrefix,ep,fh,ph(1,2:end),refEp,flip);%Second, the rest:
        
        set(ph(:,1),'CLim',[-1 1]*0.5);
        %     set(ph(:,2:end),'YTickLabels',{},'CLim',[-1 1]*2);
        set(ph(:,2:end),'YTickLabels',{},'CLim',[-1 1]*0.5);
        set(ph,'FontSize',8)
        pos=get(ph(1,end),'Position');
        axes(ph(1,end))
        colorbar
        set(ph(1,end),'Position',pos);
        set(gcf,'color','w');
        
        resDir = [scriptDir '/RegressionAnalysis/RegModelResults_V14/GroupResults/'];
        if (saveResAndFigure)
            if not(isfolder(resDir))
                mkdir(resDir)
            end
            saveas(fh, [resDir groupID '_group_Checkerboard_ver' num2str(usefft) num2str(normalizeData) regModelVersion 'split_' num2str(splitCount) 'flip_' num2str(flip) '.png'])
            %                 saveas(fh, [resDir groupID '_group_Checkerboard_ver' num2str(usefft) num2str(normalizeData) regModelVersion 'split_' num2str(splitCount) 'flip_' num2str(flip)], 'epsc')
        end
        
        if flip ~= 2 %run regression on the full (not asymmetry) data
            % run regression and save results
            format compact % format loose %(default)
            % not normalized first, then normalized, arugmnets order: (Data, normalizeData, isGroupData, dataId, resDir, saveResAndFigure, version, usefft)
            [trans1,trans2]=runRegression_Generalization(Data, false, true, [groupID regModelVersion 'split_' num2str(splitCount) 'flip_' num2str(flip)], resDir, saveResAndFigure, regModelVersion, usefft)
            runRegression_Generalization(Data, true, true, [groupID regModelVersion 'split_' num2str(splitCount) 'flip_' num2str(flip)], resDir, saveResAndFigure, regModelVersion, usefft)
        else
            asymCos = findCosBtwAsymOfEpochs(Data, size(newLabelPrefix,2))
        end
    end
    
    Rsquared1_Adj(1)=trans1.Rsquared.Adjusted;
    Rsquared1_Ord(1)=trans1.Rsquared.Ordinary;
    
    Rsquared2_Adj(1)=trans2.Rsquared.Adjusted;
    Rsquared2_Ord(1)=trans2.Rsquared.Ordinary;
    
    Coefficients_trans1(:,1)=trans1.Coefficients.Estimate;
    Coefficients_trans2(:,1)=trans2.Coefficients.Estimate;
end
% end

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
refEpOG = defineReferenceEpoch('OGbase',ep);
refEpTR = defineReferenceEpoch('TM base',ep);
post1ep = ep(strcmp(ep.Properties.ObsNames,'Post1_{Early}'),:);
post2ep = ep(strcmp(ep.Properties.ObsNames,'Post2_{Early}'),:);
post1lateep = ep(strcmp(ep.Properties.ObsNames,'Post1_{Late}'),:);
post2lateep = ep(strcmp(ep.Properties.ObsNames,'Post2_{Late}'),:);

Data = cell(1,4);

for flip = [1]%,2] %2 legs separate first (flip = 1) and then asymmetry (flip = 2)
    %     the flip asymmetry plots average of summation and the average of
    %     asymmetry.
    for i = 1:n_subjects
        if plotGroup
            adaptDataSubject = normalizedTMFullAbrupt;
            figureSaveId = groupID;
        else
            adaptDataSubject = normalizedTMFullAbrupt.adaptData{1, i};
            figureSaveId = subID{i};
        end
        
        fh=figure('Units','Normalized','OuterPosition',[0 0 1 1]);
        ph=tight_subplot(1,8)%,[.03 .005],.04,.04);
        
        % plot after effects only
        [~,~,~,Data{1},~]= adaptDataSubject.plotCheckerboards(newLabelPrefix,refEpOG,fh,ph(1,1),[],flip); %plot OG base
        [~,~,~,Data{2},~]= adaptDataSubject.plotCheckerboards(newLabelPrefix,refEpTR,fh,ph(1,2),[],flip); %plot OG base with shoe
        
        %         if (contains(groupID, 'TS'))%correct post 1 with OG no nimbus, post 2 with OG with nimbus, i.e.,TR
        %             adaptDataSubject.plotCheckerboards(newLabelPrefix,post1ep,fh,ph(1,3),refEpOG,flip); %post1 is OG
        %             adaptDataSubject.plotCheckerboards(newLabelPrefix,post2ep,fh,ph(1,4),refEpTR,flip); %post2 is with Nimbus(TR)
        %         elseif (contains(groupID, 'TR')) %correct post 1 with TR(nimbus), post 2 with OG (No nimbus)
        %             adaptDataSubject.plotCheckerboards(newLabelPrefix,post1ep,fh,ph(1,3),refEpTR,flip);
        %             adaptDataSubject.plotCheckerboards(newLabelPrefix,post2ep,fh,ph(1,4),refEpOG,flip);
        %         end
        
        [~,~,~,Data{3},~]= adaptDataSubject.plotCheckerboards(newLabelPrefix,post1ep,fh,ph(1,3),[],flip);
        [~,~,~,Data{4},~]= adaptDataSubject.plotCheckerboards(newLabelPrefix,post2ep,fh,ph(1,4),[],flip);
        [~,~,~,Data{5},~]=adaptDataSubject.plotCheckerboards(newLabelPrefix,post1lateep,fh,ph(1,5),[],flip);
        [~,~,~,Data{6},~]=adaptDataSubject.plotCheckerboards(newLabelPrefix,post2lateep,fh,ph(1,6),[],flip);
        
        %         if strcmpi(regModelVersion, 'TR')
        adaptDataSubject.plotCheckerboards(newLabelPrefix,post1lateep,fh,ph(1,7),refEpTR,flip);
        title('Post1_{Late}-TRbase_{late}')
        
        disp('Cosine(Post1_{Late},TRbase_{late})')
        cos= Cosine2Matrix(Data{5},Data{2})
        
        
        
        adaptDataSubject.plotCheckerboards(newLabelPrefix,post2lateep,fh,ph(1,8),refEpOG,flip);
        title('Post2_{Late}-OGbase_{late}')
        
        
        disp('Cosine(Post2_{Late},OGbase_{late})')
        cos= Cosine2Matrix(Data{1},Data{6})
        
        %         else
        %
        %             adaptDataSubject.plotCheckerboards(newLabelPrefix,post1lateep,fh,ph(1,7),refEpOG,flip);
        %             title('Post1_{Late}-OGbase_{late}')
        %             disp('Cosine(Post1_{Late},OGbase_{late})')
        %             cos= Cosine2Matrix(Data{1},Data{5})
        %
        %             %             cosine(Data{1},Data{3})
        %
        %             adaptDataSubject.plotCheckerboards(newLabelPrefix,post2lateep,fh,ph(1,8),refEpTR,flip);
        %             title('Post2_{Late}-TRbase_{late}')
        %             disp('Cosine(Post2_{Late},TRbase_{late})')
        %             cos= Cosine2Matrix(Data{2},Data{6})
        %             %             cosine(Data{2},Data{4})
        %         end
        
        
        %         adaptDataSubject.plotCheckerboards(newLabelPrefix,post1ep,fh,ph(1,7),ep(4,:),flip);
        %         title('trans2');
        %         adaptDataSubject.plotCheckerboards(newLabelPrefix,post2ep,fh,ph(1,8),post1lateep,flip);
        %         title('trans2');
        %         adaptDataSubject.plotCheckerboards(newLabelPrefix,refEpOG,fh,ph(1,9),refEpTR,flip);
        %         title('OG-TRbase');
        
        set(ph(:,1),'CLim',[-1 1]);
        set(ph(:,2:end),'YTickLabels',{},'CLim',[-1 1]);
        set(ph,'FontSize',8)
        pos=get(ph(1,end),'Position');
        axes(ph(1,end))
        colorbar
        set(ph(1,end),'Position',pos);
        set(gcf,'color','w');
        
        currLabelPrefix = newLabelPrefix;
        
        
        
        if (saveResAndFigure)
            resDir = [scriptDir '/RegressionAnalysis/RegModelResults_V15/'];
            if plotGroup
                resDir = [resDir 'GroupResults/'];
            end
            if not(isfolder(resDir))
                mkdir(resDir)
            end
            if flip == 1
                saveas(fh, [resDir figureSaveId '_CheckerBoard_AE_NoBase.png'])
                %                 saveas(fh, [resDir figureSaveId '_AllEpochCheckerBoard_' num2str(session)],'epsc')
            else
                saveas(fh, [resDir figureSaveId '_CheckerBoard_AE_NoBase_Asym.png'])
            end
        end
        
        if plotGroup
            break
        end
    end
end



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