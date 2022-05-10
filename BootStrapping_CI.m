%% Bootstraping data 
%% Load data and Plot checkerboard for all conditions.
clear; close all; clc;

% set script parameters, SHOULD CHANGE/CHECK THIS EVERY TIME.
groupID = 'PATS';
saveResAndFigure = false;
plotAllEpoch = true;
plotIndSubjects = true;
plotGroup = true;
bootstrap=true;

scriptDir = fileparts(matlab.desktop.editor.getActiveFilename);
files = dir ([scriptDir '/data/' groupID '*params.mat']);

% n_subjects = size(files,1);
n_subjects = size(files,1);
% p = randperm(4,4);
ii=0;

subID = cell(1, n_subjects);
sub=cell(1,n_subjects);
% for i = 1:n_subjects
%     sub{i} = files(i).name;
%     subID{i} = sub{i}(1:end-10);
%    
% end

for i =1:n_subjects
    ii=1+ii;
    sub{ii} = files(i).name;
    subID{ii} = sub{ii}(1:end-10);
end

subID

regModelVersion =  'default';
%% load and prep data

muscleOrder={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT', 'VM', 'VL', 'RF', 'TFL', 'HIP', 'GLU'};
n_muscles = length(muscleOrder);

ep=defineEpochs_regressionYA('nanmean');
refEpTM = defineReferenceEpoch('TM base',ep);
refEpOG = defineReferenceEpoch('OG base',ep);

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
Adaptearly=defineReferenceEpoch('Adaptation_{early}',ep);
Pos1_Late=defineReferenceEpoch('Post1_{Late}',ep);
Pos2_Early=defineReferenceEpoch('Post2_{Early}',ep);

if (contains(groupID, 'ATR'))
    AfterPos= defineReferenceEpoch('TM_{afterPos}',ep);
    AfterNeg=defineReferenceEpoch('TM_{afterNeg}',ep);
    
elseif (contains(groupID, 'ATS'))
    AfterPos= defineReferenceEpoch('OG_{afterPos}',ep);
    AfterNeg=defineReferenceEpoch('OG_{afterNeg}',ep);
end

refEpPost1Early= defineReferenceEpoch('Post1_{Early}',ep);
if contains(groupID,'ATR')
    refEp= defineReferenceEpoch('TM base',ep);
elseif contains(groupID,'ATS')
    refEp= defineReferenceEpoch('OG base',ep);
end

%% Bootstrapping section
saveResAndFigure=0;
n_subjects = length(subID);

Trans1=[];
Trans2=[];
Trans3=[];


GroupData=adaptationData.createGroupAdaptData(sub); %loading the data
GroupData=GroupData.removeBadStrides; %Removing bad strides

newLabelPrefix = defineMuscleList(muscleOrder);
normalizedGroupData = GroupData.normalizeToBaselineEpoch(newLabelPrefix,refEpOG); %Normalized by OG base same as nimbus data
ll=normalizedGroupData.adaptData{1}.data.getLabelsThatMatch('^Norm');
l2=regexprep(regexprep(ll,'^Norm',''),'_s','');
normalizedGroupData=normalizedGroupData.renameParams(ll,l2);
newLabelPrefix = regexprep(newLabelPrefix,'_s','');
flip=1;
Data=cell(n_subjects,7);
summFlag='nanmedian';
% [data,~,~,groupData] = normalizedGroupData.getCheckerboardsData(newLabelPrefix,PosShort,TMbeforePos,2,summFlag);
% normalizedGroupData.plotCheckerboards(newLabelPrefix,AdaptLate,[],[],[],2,summFlag);
% normalizedGroupData.plotCheckerboards(newLabelPrefix,Adaptearly,[],[],[],2,summFlag);

 
  for s=1:n_subjects

    [Data(s,1:6),regressorNames]=getRegressionsData(newLabelPrefix,normalizedGroupData,[],1,0,NegShort,TMbeforeNeg,PosShort,...
    TMbeforePos,AdaptLate,Pos1_Early,Pos1_Late,Pos2_Early, AfterPos, OGbase, TMbase,s,flip);
    adaptDataSubject = normalizedGroupData.adaptData{1, s};
    
    Data{s,7} = adaptDataSubject.getCheckerboardsData(newLabelPrefix,refEpPost1Early,refEp,flip); %|EMG_earlyPost1 -  EMG_Baseline|
    
  end

%% Shuffeling the data - per muscle
  for r=1:7
      p = randperm(28);
      disp(['r=' num2str(r)])
      for s=1:4
                 
          Data{s,r}= Data{s,r}(:,p);
          
          
      end
  end

%%

if bootstrap
    n=2000; %number of iterations
   f = waitbar(0,'1','Name','Boostrapping Data',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

setappdata(f,'canceling',0);

for i=1:n
    
    temp=[];
    temp2={};
    
    % Check for clicked Cancel button
    if getappdata(f,'canceling')
        break
    end
    % Update waitbar and message
    ww=waitbar(i/n,f,['Iteration ' num2str(i)]);
    
    %         normalizedGroupData =[];
    DataBoot=datasample(Data,length(subID),'Replace',true);
    
    for c=1:7
        for s=1:length(subID)
            temp(:,:,s)=DataBoot{s,c};
        end
        temp2{c}=nanmedian(temp,3);
    end
    
    
    
    nw=datestr(now,'yy-mm-dd');
    resDir = [scriptDir '/RegressionAnalysis/RegModelResults_',nw ,'/GroupResults/'];
    
    
    % run regression and save results
    format compact % format loose %(default)
    % not normalized first, then normalized, arugmnets order: (Data, normalizeData, isGroupData, dataId, resDir, saveResAndFigure, version, usefft)
    [trans1,trans2,trans3, vec_norm] = runRegression_Generalization(temp2, false, true,...
        [groupID regModelVersion 'split_' num2str(splitCount) 'flip_' num2str(flip)], resDir,...
        saveResAndFigure, regModelVersion, usefft);
    %            [trans1Norm,trans2Norm,tran3Norm] = runRegression_Generalization(Data, true, true, [groupID regModelVersion 'split_' num2str(splitCount) 'flip_' num2str(flip)], resDir, saveResAndFigure, regModelVersion, usefft);
    
    %Trans 1
    Trans1.Rsquared_Ord(i,1)= trans1.Rsquared.Ordinary;
    Trans1.Rsquared_Adj(i,1) = trans1.Rsquared.Adjusted;
    Trans1.betas(i,:) = trans1.Coefficients.Estimate';
    Norm_Regressors(i,:)= vec_norm;
    
    %Trans 2
    Trans2.Rsquared_Ord(i,1)= trans2.Rsquared.Ordinary;
    Trans2.Rsquared_Adj(i,1) = trans2.Rsquared.Adjusted;
    Trans2.betas(i,:) = trans2.Coefficients.Estimate';
    
    %Trans 2
    Trans3.Rsquared_Ord(i,1)= trans3.Rsquared.Ordinary;
    Trans3.Rsquared_Adj(i,1) = trans3.Rsquared.Adjusted;
    Trans3.betas(i,:) = trans3.Coefficients.Estimate';
    
    
    
%     summFlag='nanmedian';
%     [Data2] = normalizedGroupData.getCheckerboardsData(newLabelPrefix,refEpPost1Early,refEp,flip,summFlag); %|EMG_earlyPost1 -  EMG_Baseline|
%     Data2 = nanmedian(Data2, 4);
    Norm_AF(i)= norm(reshape(temp2{7},[],1));
    
    
    
   
end
     close(ww)
    
end

delete(f)

% nw='22-01-20';
resDir = [scriptDir '/RegressionAnalysis/RegModelResults_', nw,'/BootstrappingResults/'];
 
if not(isfolder(resDir))
    mkdir(resDir)
end
% save([resDir groupID '_group_iterations_' num2str(i) '_numberOfSub_' num2str(n_subjects)], 'Trans1','Trans2','Trans3','Norm_Regressors','Norm_AF')
save([resDir groupID '_ShuffleDataBySubjectV3_group_iterations_' num2str(i) '_numberOfSub_' num2str(n_subjects)], 'Trans1','Trans2','Trans3','Norm_Regressors','Norm_AF','subID')