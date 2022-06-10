%Getting plot of the epochs of interest this will help us to see the
%quality of the data 

%% Load data and Plot checkerboard for all conditions.
clear; close all; clc;

% set group/participant ID that you want to plot 
groupID = 'PATS';

% scriptDir = fileparts(matlab.desktop.editor.getActiveFilename %you can
% use the path format but I am currently working on the folder - to avoid
% colling other versions of the params files
files = dir ([ groupID '*params.mat']);

n_subjects = size(files,1);


subID = cell(1, n_subjects);
sub=cell(1,n_subjects);

ii=0;
for i =1:n_subjects
    ii=1+ii;
    sub{ii} = files(i).name;
    subID{ii} = sub{ii}(1:end-10);
end

subID


% Defining the muscle and the order they are going to be plot 
muscleOrder={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT', 'VM', 'VL', 'RF','TFL', 'GLU'};
n_muscles = length(muscleOrder);

%Defining epochs of interest 
ep=defineROI_Exploration('nanmean');
%define epoch to normalize the data against
refEpTM = defineReferenceEpoch('TMbase',ep);

%% EMG normalization 
GroupData=adaptationData.createGroupAdaptData(sub); %loading the data
GroupData=GroupData.removeBadStrides; %Removing bad strides


newLabelPrefix = defineMuscleList(muscleOrder); %name of the muscles 
normalizedGroupData = GroupData.normalizeToBaselineEpoch(newLabelPrefix,refEpTM); %Normalized by OG base same as nimbus data

%Renaming the muscles labels - we remove the normalization label
ll=normalizedGroupData.adaptData{1}.data.getLabelsThatMatch('^Norm'); 
l2=regexprep(regexprep(ll,'^Norm',''),'_s','');
normalizedGroupData=normalizedGroupData.renameParams(ll,l2);
newLabelPrefix = regexprep(newLabelPrefix,'_s','');


%% Plotting all epochs of interest 

summFlag='nanmedian'; %how are we going to get the group behavior, due to the high noise on the data we choose median instead of mean
fs=8; %font size, this is for the plots 
flip=2; %1 individual leg behavior - 2: Asymmetry value 

if flip==1
    n=2;
    method='IndvLegs';
else
    n=1;
    method='Asym';
    % Color definition - This color map will match the color map form the
    % simulaiton plots
    ex1=[1,0,0];
    ex2=[0,0,1];
    cc=[0    0.4470    0.7410
        0.8500    0.3250    0.0980
        0.9290    0.6940    0.1250
        0.4940    0.1840    0.5560
        0.4660    0.6740    0.1880
        0.3010    0.7450    0.9330
        0.6350    0.0780    0.1840];
    ex1=cc(2,:);
    ex2=cc(5,:);
    mid=ones(1,3);
    N=100;
    gamma=1.5; %gamma > 1 expands the white (mid) part of the map, 'hiding' low values. Gamma<1 does the opposite
    gamma=1;
    map=[flipud(mid+ (ex1-mid).*([1:N]'/N).^gamma); mid; (mid+ (ex2-mid).*([1:N]'/N).^gamma)];
    
    
end

%plotting epochs of interest bias 
fh=figure('Units','Normalized','OuterPosition',[0 0 1 1]);
ph=tight_subplot(1,length(ep),[.03 .005],.04,.04);
normalizedGroupData.plotCheckerboards(newLabelPrefix,ep,fh,ph,[],flip,summFlag);
set(ph(:,1),'CLim',[-1 1]*1,'FontSize',fs);
set(ph(:,2:end),'YTickLabels',{},'CLim',[-1 1]*1,'FontSize',fs);

%plotting epochs of interest unbias - we are removing TM base for all
%condition - NO smart bc we have OG data too - TO BE FIX 
fh=figure('Units','Normalized','OuterPosition',[0 0 1 1]);
ph=tight_subplot(1,length(ep),[.03 .005],.04,.04);
normalizedGroupData.plotCheckerboards(newLabelPrefix,ep,fh,ph,refEpTM,flip,summFlag);
set(ph(:,1),'CLim',[-1 1]*1,'FontSize',fs);
set(ph(:,2:end),'YTickLabels',{},'CLim',[-1 1]*1,'FontSize',fs);
%% Plotting eopochs of interest next to each other to be able to visually compare between them
unbias=0;
if unbias==1
    bias = [];
else
    bias = defineReferenceEpoch('TMbase',ep);
end
cc=[];

%Positive hort perturbations behavior 
cc{1}={'PostShort_{early}','PosShort_{late}','PosRamp','Adaptation_{early}'};
%Baseline 
cc{2}={'TMfast','TMslow','TMmid1','TMbase','OGbase','Post1_{late}','Post2_{late}','TMmid2_{early}','OG2_{early}'};
%Epoch of interest for the models - This can give us inside in which muscle
%are the problematics on the fitting 
cc{3}={'OGbase','TMbase','Adaptation_{early}','Adaptation_{mid}','Adaptation_{late}','Post1_{early}','Optimal'};
%This is mostly to evalute the post adaptation  - we are having issue
%reconstruction this data.
cc{4}={'PostShort_{early}','PosShort_{late}','PosRamp','NegShort_{early}','NegShort_{late}','Post1_{early}'};


for c=1:length(cc)
fh=figure('Units','Normalized','OuterPosition',[0 0 1 1]);
ph=tight_subplot(1,length(cc{c}),[.03 .005],.04,.04);

for i=1:length(cc{c})
    ep2=defineReferenceEpoch(cc{c}{i},ep);
    normalizedGroupData.plotCheckerboards(newLabelPrefix,ep2,fh,ph(1,i),bias,flip,summFlag);
end
set(ph(:,1),'CLim',[-1 1]*1,'FontSize',fs);
set(ph(:,2:end),'YTickLabels',{},'CLim',[-1 1]*1,'FontSize',fs);
if flip==2
colormap(flipud(map))
end
colorbar   
end

