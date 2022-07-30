%getting Cs 
%% Load data and Plot checkerboard for all conditions.
clear; close all; clc;

% set script parameters, SHOULD CHANGE/CHECK THIS EVERY TIME.
groupID = 'PATR';

% scriptDir = fileparts(matlab.desktop.editor.getActiveFilename);
files = dir ([ groupID '*params.mat']);


n_subjects = size(files,1);

ii=0;
for i =1:n_subjects
    ii=1+ii;
    sub{ii} = files(i).name;
    subID{ii} = sub{ii}(1:end-10);
end

subID

regModelVersion =  'default';

%%%% load and prep data
muscleOrder={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT', 'VM', 'VL', 'RF', 'HIP','TFL', 'GLU'};

n_muscles = length(muscleOrder);


ep=defineRegressorsDynamics('nanmean');
refEpTM = defineReferenceEpoch('TM base',ep);


%%

GroupData=adaptationData.createGroupAdaptData(sub); %loading the data
GroupData=GroupData.removeBadStrides; %Removing bad strides

% GroupData=group;

newLabelPrefix = defineMuscleList(muscleOrder);
normalizedGroupData = GroupData.normalizeToBaselineEpoch(newLabelPrefix,refEpTM); %Normalized by OG base same as nimbus data
ll=normalizedGroupData.adaptData{1}.data.getLabelsThatMatch('^Norm');
l2=regexprep(regexprep(ll,'^Norm',''),'_s','s');
normalizedGroupData=normalizedGroupData.renameParams(ll,l2);
newLabelPrefix = regexprep(newLabelPrefix,'_s','s');

Data=cell(7,1);
Data2=cell(5,1);
group=cell(5,1);
summFlag='nanmedian';

%% Remove aftereffects using Shuqi's code

 removeBadmuscles=0;
if removeBadmuscles==1
     [RemovedData]=RemoveBadMuscles(normalizedGroupData,groupID);
    normalizedGroupData=RemovedData;
end
%%

epochOfInterest={'TM base','TM mid 1','PosShort_{early}','PosShort_{late}','Ramp','Optimal'};
fh=figure('Units','Normalized','OuterPosition',[0 0 1 1]);
ph=tight_subplot(1,length(epochOfInterest),[.03 .005],.04,.04);

flip=1;

if flip==1
n=2;
method='IndvLegs';
else
   n=1; 
   method='Asym';
end

C=[];
for l=1:length(epochOfInterest)
ep2=defineReferenceEpoch(epochOfInterest{l},ep);
normalizedGroupData.plotCheckerboards(newLabelPrefix,ep2,fh,ph(1,l),[],flip,summFlag);
[~,~,~,Data{l}]=normalizedGroupData.getCheckerboardsData(newLabelPrefix,ep2,[],flip,summFlag);
C=[C reshape(Data{l}(:,end:-1:1),12*n_muscles*n,1)];
end 
  
%%
% resDir = [cd '/LTI models/'];
% save([resDir '/'  groupID,'_',num2str(n_subjects),'_',method,'C',num2str(length(epochOfInterest)) ,'_ShortPertubations_RemovedBadMuscle_',num2str(removeBadmuscles)], 'C','epochOfInterest')

%%
%% Color definition 
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

%%
% colormap(map)
fs=14;
colormap(flipud(map))
% colormap default
set(gcf,'color','w');
colorbar                                                                                                                                                                                         
set(ph(:,1),'CLim',[-1 1]*1,'FontSize',fs);
set(ph(:,2:end),'YTickLabels',{},'CLim',[-1 1]*1,'FontSize',fs);

%%
