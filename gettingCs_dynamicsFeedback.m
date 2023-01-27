%getting Cs
%% Load data and Plot checkerboard for all conditions.
clear; close all; clc;

groupID ='BAT';
[normalizedGroupData, newLabelPrefix,n_subjects]=creatingGroupdataWnormalizedEMG(groupID);

Data=cell(7,1);
Data2=cell(5,1);
group=cell(5,1);
summFlag='nanmedian';
n_muscles =14; %number fo muscle in the experiment
%% Remove aftereffects

removeBadmuscles=0;
if removeBadmuscles==1

    normalizedGroupData= RemovingBadMuscleToSubj(normalizedGroupData);
end
%% Getting the data that we use for matrices
ep=defineRegressorsDynamicsFeedback('nanmean'); % This is a hardcade (aka specific to this experiment) script to get the data from name epochs of interest

epochOfInterest={'TM base','TM mid 1','PosShort_{early}',...
    'PosShort_{late}','Ramp','Optimal','Adaptation',...
    'Adaptation_{early}','TiedPostPos','OG2','NegShort_{late}',...
    'Post1_{Early}','TMbase_{early}'}; % This line chooses the epochs that want to get data from


% epochOfInterest={'Post1_{Early}','TiedPostPos','TMmid2','Tied post ramp','Ramp','Tied post Neg'};
% epochOfInterest={'OG2','Ramp','Optimal'};


fh=figure('Units','Normalized','OuterPosition',[0 0 1 1]);
ph=tight_subplot(1,length(epochOfInterest),[.03 .005],.04,.04);

flip=2; %1 for individual leg analysis and 2 for asymmetric (R-L)
if flip==1 
    n=2;
    method='IndvLegs';
else
    n=1;
    method='Asym';
end
fdr=.1;
C=[];
for l=1:length(epochOfInterest)
    ep2=defineReferenceEpoch(epochOfInterest{l},ep); 
    [~,~,~,Data2{l}]=normalizedGroupData.plotCheckerboards(newLabelPrefix,ep2,fh,ph(1,l),[],flip,summFlag); %plotting the data
    [~,~,~,Data{l}]=normalizedGroupData.getCheckerboardsData(newLabelPrefix,ep2,[],flip,summFlag); %getting the data from the plots
    C=[C reshape(Data{l}(:,end:-1:1),12*n_muscles*n,1)]; %Reshaping the data form the plot to a vector
end

%% Saving the data
resDir = [cd];% '/LTI models/'];
save([resDir '/'  groupID,'_',num2str(n_subjects), '_',method,'C',num2str(length(epochOfInterest)) ,'_ShortPertubations_RemovedBadMuscle_',num2str(removeBadmuscles)], 'C','epochOfInterest')

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
    0.6350      0.0780    0.1840];
ex1=cc(2,:);
ex2=cc(5,:);
mid=ones(1,3);
N=100;
gamma=1.5; %gamma > 1 expands the white (mid) part of the map, 'hiding' low values. Gamma<1 does the opposite
gamma=1;
map=[flipud(mid+ (ex1-mid).*([1:N]'/N).^gamma); mid; (mid+ (ex2-mid).*([1:N]'/N).^gamma)];

%% Making the figure a bit prettier 
fs=14; %font size 
colormap(flipud(map)) %changing the color map to the one tha defined about, we flipup the matrix bc the code does L-R and we want R-L 
set(gcf,'color','w'); %setting the background white 
set(ph(:,1),'CLim',[-1 1]*1,'FontSize',fs); %making sure that the first plot color scheme goes [-1 1] and making the name of the labels larger 
set(ph(:,2:end),'YTickLabels',{},'CLim',[-1 1]*1,'FontSize',fs);
colorbar %Showing the colormap bar 
