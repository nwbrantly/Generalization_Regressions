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
removeBadmuscles=1;
if removeBadmuscles==1
    
    normalizedGroupData= RemovingBadMuscleToSubj(normalizedGroupData);
end

%% Pick muscles that you wanted to get the data from

muscleOrder={'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT'};
newLabelPrefix2 = defineMuscleList(muscleOrder); %List of muscle
newLabelPrefix2 = regexprep(newLabelPrefix2,'_s','s'); %removing the underscore "_"

for m=1:length(newLabelPrefix2)
    
    idx(m)=find(strcmp(newLabelPrefix,newLabelPrefix2(m))); %getting the index of the muscles
    
end
wanted_Muscles= newLabelPrefix(sort(idx)); %It needs to be the muscle form the slow leg first
newLabelPrefix= wanted_Muscles;
n_muscles =length(muscleOrder); %number fo muscle in the experiment
%% Getting the data that we use for matrices

ep=defineRegressorsDynamicsFeedback('nanmean'); % This is a hardcade (aka specific to this experiment) script to get the data from name epochs of interest
epochOfInterest={'TM base','TM mid 1','PosShort_{early}',...
    'PosShort_{late}','Ramp','Optimal','Adaptation',...
    'Adaptation_{early}','TiedPostPos','OG2','NegShort_{early}','NegShort_{late}',...
    'Post1_{Early}','TMbase_{early}','Tied post Neg','Tied post ramp','OG base'}; % This line chooses the epochs that want to get data from


% epochOfInterest={'Post1_{Early}','TiedPostPos','TMmid2','Tied post ramp','Ramp','Tied post Neg'};
% epochOfInterest={'Ramp','Optimal'};
% epochOfInterest={'TM base','Ramp','Optimal','Tied post ramp'}; % This line chooses the epochs that want to get data from

fh=figure('Units','Normalized','OuterPosition',[0 0 1 1]);
ph=tight_subplot(1,length(epochOfInterest),[.03 .005],.04,.04);

flip=1; %1 for individual leg analysis and 2 for asymmetric (R-L)
if flip==1
    n=2;
    method='IndvLegs';
else
    n=1;
    method='Asym';
end
fdr=.1;
C=[];

removeBias=0; %Flag for bias removal 
base=defineReferenceEpoch('TM base',ep); %pick the condition for baseline 

for l=1:length(epochOfInterest)
    ep2=defineReferenceEpoch(epochOfInterest{l},ep);
    
    if removeBias==1
        [~,~,~,Data2{l}]=normalizedGroupData.plotCheckerboards(newLabelPrefix,ep2,fh,ph(1,l),base,flip,summFlag); %plotting the data
        [~,~,~,Data{l}]=normalizedGroupData.getCheckerboardsData(newLabelPrefix,ep2,base,flip,summFlag); %getting the data from the plots
    else
        [~,~,~,Data2{l}]=normalizedGroupData.plotCheckerboards(newLabelPrefix,ep2,fh,ph(1,l),[],flip,summFlag); %plotting the data
        [~,~,~,Data{l}]=normalizedGroupData.getCheckerboardsData(newLabelPrefix,ep2,[],flip,summFlag); %getting the data from the plots
    end
    
    C=[C reshape(Data{l}(:,end:-1:1),12*n_muscles*n,1)]; %Reshaping the data form the plot to a vNoah tor
    
end

%% Saving the data
resDir = [cd];% '/LTI models/'];
save([resDir '/'  groupID,'_',num2str(n_subjects), '_',method,'C',num2str(length(epochOfInterest)) ,'_ShortPertubations_RemovedBadMuscle_',num2str(removeBadmuscles), 'RemoveBias_',num2str(removeBias),'_PosteriorMuscles'], 'C','epochOfInterest')

%%
%% Color definition
ex2=[0.2314    0.2980    0.7529];
ex1=[0.7255    0.0863    0.1608];
% ex1=[1,0,0];
% ex2=[0,0,1];
cc=[0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350      0.0780    0.1840];

mid=ones(1,3);
N=100;
% gamma=1.5; %gamma > 1 expands the white (mid) part of the map, 'hiding' low values. Gamma<1 does the opposite
gamma=1;
map=[flipud(mid+ (ex1-mid).*([1:N]'/N).^gamma); mid; (mid+ (ex2-mid).*([1:N]'/N).^gamma)];

%% Making the figure a bit prettier
fs=14; %font size
colormap(flipud(map)) %changing the color map to the one tha defined about, we flipup the matrix bc the code does L-R and we want R-L
% c=flipud('gray');
% colormap(flipud(gray))
set(gcf,'color','w'); %setting the background white
set(ph(:,1),'CLim',[-1 1]*1,'FontSize',fs); %making sure that the first plot color scheme goes [-1 1] and making the name of the labels larger
set(ph(:,2:end),'YTickLabels',{},'CLim',[-1 1]*1,'FontSize',fs);
colorbar %Showing the colormap bar
