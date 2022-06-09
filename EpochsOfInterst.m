%Getting plot of the epochs of interest this will help us to see the
%quality of the data 

%% Load data and Plot checkerboard for all conditions.
clear; close all; clc;

% set script parameters, SHOULD CHANGE/CHECK THIS EVERY TIME.
groupID = 'PATR03';
saveResAndFigure = false;
plotAllEpoch = true;
plotIndSubjects = true;
plotGroup = true;
bootstrap=true;

% scriptDir = fileparts(matlab.desktop.editor.getActiveFilename);
files = dir ([ groupID '*params.mat']);

n_subjects = size(files,1);

ii=0;

subID = cell(1, n_subjects);
sub=cell(1,n_subjects);


for i =1:n_subjects
    ii=1+ii;
    sub{ii} = files(i).name;
    subID{ii} = sub{ii}(1:end-10);
end

subID

regModelVersion =  'default';

%%%% getting data
muscleOrder={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT', 'VM', 'VL', 'RF','TFL', 'GLU'};
n_muscles = length(muscleOrder);
ep=defineROI_Exploration('nanmean');
refEpTM = defineReferenceEpoch('TMbase',ep);

%%

GroupData=adaptationData.createGroupAdaptData(sub); %loading the data
GroupData=GroupData.removeBadStrides; %Removing bad strides


newLabelPrefix = defineMuscleList(muscleOrder);
normalizedGroupData = GroupData.normalizeToBaselineEpoch(newLabelPrefix,refEpTM); %Normalized by OG base same as nimbus data
ll=normalizedGroupData.adaptData{1}.data.getLabelsThatMatch('^Norm');
l2=regexprep(regexprep(ll,'^Norm',''),'_s','');
normalizedGroupData=normalizedGroupData.renameParams(ll,l2);
newLabelPrefix = regexprep(newLabelPrefix,'_s','');


summFlag='nanmedian';

%%
epochOfInterest={'Adaptation_{early}','Adaptation','Post1_{Early}','NegShort_{early}','TM base', };
% fh=figure('Units','Normalized','OuterPosition',[0 0 1 1]);
% ph=tight_subplot(1,length(ep),[.03 .005],.04,.04);

fs=10;
flip=1;
if flip==1
n=2;
method='IndvLegs';
else
   n=1; 
   method='Asym';
end

fh=figure('Units','Normalized','OuterPosition',[0 0 1 1]);
ph=tight_subplot(1,length(ep),[.03 .005],.04,.04);
normalizedGroupData.plotCheckerboards(newLabelPrefix,ep,fh,ph,[],flip,summFlag);
set(ph(:,1:end),'YTickLabels',{},'CLim',[-1 1]*1,'FontSize',fs);
fh=figure('Units','Normalized','OuterPosition',[0 0 1 1]);
ph=tight_subplot(1,length(ep),[.03 .005],.04,.04);
normalizedGroupData.plotCheckerboards(newLabelPrefix,ep,fh,ph,refEpTM,flip,summFlag);
set(ph(:,1:end),'YTickLabels',{},'CLim',[-1 1]*1,'FontSize',fs);

  
%%
% resDir = [cd '/LTI models/']
% resDir = [cd];
% save([resDir '/'  groupID,'_',num2str(n_subjects),'_',method,'C_WO_HIP',num2str(length(epochOfInterest))], 'C', 'epochOfInterest')

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
% fs=14;
% colormap(flipud(map))
% % colormap default
% set(gcf,'color','w');
% colorbar                                                                                                                                                                                         
% set(ph(:,1),'CLim',[-1 1]*2,'FontSize',fs);
set(ph(:,1:end),'YTickLabels',{},'CLim',[-1 1]*2,'FontSize',fs);

%%

