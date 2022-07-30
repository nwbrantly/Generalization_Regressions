%getting Cs 
%% Load data and Plot checkerboard for all conditions.
% clear; close all; clc;

% set script parameters, SHOULD CHANGE/CHECK THIS EVERY TIME.
groupID = 'ATR';
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

%%%% load and prep data

% muscleOrder={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT', 'VM', 'VL', 'RF', 'TFL', 'HIP', 'GLU'};
% muscleOrder={'GLU','TFL','HIP','RF','VL','VM','BF', 'SEMB','SEMT','MG','LG','SOL','PER','TA'};
muscleOrder={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT', 'VM', 'VL', 'RF', 'HIP','TFL', 'GLU'};
% muscleOrder={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT', 'VM', 'VL', 'RF','TFL', 'GLU'};

% muscleOrder={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT', 'VM', 'VL', 'RF', 'HIP', 'ADM', 'TFL', 'GLU'};
% muscleOrder(end:-1:1) = muscleOrder(:);

% muscleOrder={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT', 'VM', 'VL', 'RF', 'HIP', 'ADM', 'TFL', 'GLU'};
% muscleOrder={'GLU','HIP','TFL','RF','VL','VM','SEMT','SEMB','BF','MG','LG','SOL','PER','TA'};

n_muscles = length(muscleOrder);
ep=defineEpochs_regressionYA('nanmean');
refEpTM = defineReferenceEpoch('TM base',ep);

%%

GroupData=adaptationData.createGroupAdaptData(sub); %loading the data
GroupData=GroupData.removeBadStrides; %Removing bad strides

% GroupData=group;

newLabelPrefix = defineMuscleList(muscleOrder);
normalizedGroupData = GroupData.normalizeToBaselineEpoch(newLabelPrefix,refEpTM); %Normalized by OG base same as nimbus data
ll=normalizedGroupData.adaptData{1}.data.getLabelsThatMatch('^Norm');
l2=regexprep(regexprep(ll,'^Norm',''),'_s','');
normalizedGroupData=normalizedGroupData.renameParams(ll,l2);
newLabelPrefix = regexprep(newLabelPrefix,'_s','');

Data=cell(7,1);
Data2=cell(5,1);
group=cell(5,1);
summFlag='nanmedian';

% epochOfInterest={'TM base','Adaptation_{early}','Adaptation','Post1_{Early}','NegShort_{early}'};
% epochOfInterest={'Post1_{Early}','NegShort_{early}'};
%%
epochOfInterest={'Adaptation_{early}','Adaptation','Post1_{Early}','NegShort_{early}','TM base', };
fh=figure('Units','Normalized','OuterPosition',[0 0 1 1]);
ph=tight_subplot(1,length(epochOfInterest),[.03 .005],.04,.04);

flip=2;
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
% [~,~,~,Data{l}]=
normalizedGroupData.plotCheckerboards(newLabelPrefix,ep2,fh,ph(1,l),[],flip,summFlag);
[ref,~,~,Data{l}]=normalizedGroupData.getCheckerboardsData(newLabelPrefix,ep2,[],flip,summFlag);

% [dataPoint]=getEarlyLateData_v2(this,labels,ep2.Condition,removeBiasFlag,numberOfStrides(i),exemptLast(i),exemptFirst(i),padWithNaNFlag);%Get data
    C=[C reshape(Data{l}(:,end:-1:1),12*n_muscles*n,1)];
end 





% AdaptationEarly=Data{1}(:,end:-1:1);
% C1=reshape(AdaptationEarly,12*n_muscles*n,1);
% AdaptationLate=Data{2}(:,end:-1:1);
% C2=reshape(AdaptationLate,12*n_muscles*n,1);
% if length(epochOfInterest)==3
%     NegShort=Data{3}(:,end:-1:1);
%     C3=reshape(NegShort,12*n_muscles*n,1);
%     C=[C1 C2 C3];
%     
% elseif length(epochOfInterest)==4
%     NegShort=Data{3}(:,end:-1:1);
%     C3=reshape(NegShort,12*n_muscles*n,1);
%     Post1=Data{4}(:,end:-1:1);
%     C4=reshape(Post1,12*n_muscles*n,1);
%     C=[C1 C2 C3 C4];
% elseif length(epochOfInterest)==5
%     NegShort=Data{3}(:,end:-1:1);
%     C3=reshape(NegShort,12*n_muscles*n,1);
%     Post1=Data{4}(:,end:-1:1);
%     C4=reshape(Post1,12*n_muscles*n,1);
%     TMbase=Data{5}(:,end:-1:1);
%     C5=reshape(TMbase,12*n_muscles*n,1);
%     C=[C1 C2 C3 C4 C5];
% else
%     C=[C1 C2];
% end

  
%%
resDir = [cd '/LTI models_ATS_ATR/']
% resDir = [cd];
save([resDir '/'  groupID,'_',num2str(n_subjects),'_',method,'C',num2str(length(epochOfInterest))], 'C', 'epochOfInterest')

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
set(ph(:,1),'CLim',[-1 1]*2,'FontSize',fs);
set(ph(:,2:end),'YTickLabels',{},'CLim',[-1 1]*2,'FontSize',fs);

%%
ph=figure 


% ytl={'GLU','TFL','HIP','RF','VL','VM','BF', 'SEMB','SEMT','MG','LG','SOL','PER','TA'};
ytl={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT', 'VM', 'VL', 'RF', 'HIP','TFL', 'GLU'};
% ytl=newLabelPrefix;
ytl=ytl(end:-1:1);
yt=1:14*2;
fs=14;

subplot(1,4,1)
imagesc((reshape(C(:,1),12,14)'))
caxis([-1 1])
set(gca,'XTick',[],'YTick',yt,'YTickLabel',ytl,'FontSize',fs)
title('Adaptation_{early}')
% title('F(A)')


subplot(1,4,2)
imagesc((reshape(C(:,2),12,14)'))
caxis([-1 1])
set(gca,'XTick',[],'YTick',yt,'YTickLabel',ytl,'FontSize',fs)
title('Adaptation_{late}')
% title(['F(-A)'])

subplot(1,4,3)
imagesc((reshape(C(:,3),12,14)'))
caxis([-1 1])
set(gca,'XTick',[],'YTick',yt,'YTickLabel',ytl,'FontSize',fs)
% title(['Negative_{Short}'])

%%

           %Colormap:
%             ex2=[0.2314    0.2980    0.7529];
%             ex1=[0.7255    0.0863    0.1608];
%             gamma=.5;
%             map=[bsxfun(@plus,ex1.^(1/gamma),bsxfun(@times,1-ex1.^(1/gamma),[0:.01:1]'));bsxfun(@plus,ex2.^(1/gamma),bsxfun(@times,1-ex2.^(1/gamma),[1:-.01:0]'))].^gamma;

%             colormap(flipud(map))
%%
subplot(1,4,4)
imagesc((reshape(C(:,1),12,14)'))
caxis([-1 1])
set(gca,'XTick',[],'YTick',yt,'YTickLabel',ytl,'FontSize',fs)
title(['-(Adaptation_ /{early})'])
%%
% CosTheta = max(min(dot(C(:,1),C(:,2))/(norm(C(:,1))*norm(C(:,2))),1),-1);
% ThetaInDegrees = real(acosd(CosTheta));
% title(['Adaptation_{late}'],[ 'cosine = ' num2str(ThetaInDegrees)])

subplot(1,4,3)
imagesc((reshape(C(:,3),12,14)'))
caxis([-1 1])
set(gca,'XTick',[],'YTick',yt,'YTickLabel',ytl,'FontSize',fs)

%%
                                                                                                                                                                                      
% set(ph(:,1),'CLim',[-1 1]*1);
CosTheta = max(min(dot(C(:,1),C(:,3))/(norm(C(:,1))*norm(C(:,3))),1),-1);
ThetaInDegrees = real(acosd(CosTheta));

title(['NegShort_{early}'],[ 'cosine = ' num2str(ThetaInDegrees)])

subplot(1,4,4)
imagesc((reshape(abs(C(:,1))-abs(C(:,3)),12,14)'))
caxis([-1 1])

title('|Adaptation_{early}| - |NegShrot_{early}|')
set(gca,'XTick',[],'YTick',yt,'YTickLabel',ytl,'FontSize',fs)

colormap(flipud(map))
set(gcf,'color','w');
colorbar   

vif(C)

%% Getting data - Median across participants and then mean across strides
dataPoints=[];
Yasym=[];
dataAll={};
epochOfInterest={'Adaptation_{early}','Adaptation','Post1_{Early}','TM base'};
padWithNaNFlag=true;
removeBiasFlag=0;

[~,labels,~]=normalizedGroupData.adaptData{1}.getPrefixedEpochData(newLabelPrefix, ep(1,:),true);
labels=reshape(labels,size(labels,1)*size(labels,2),1);
labels=[labels(length(labels)/2+1:end); labels(1:length(labels)/2)];
 Yasym=[];
 Regressors=[];
%%
for l=1:length(epochOfInterest)
    ep2=defineReferenceEpoch(epochOfInterest{l},ep);
    
     aux1=[];
    %Second: get data
    conds=ep2.Condition;
    numberOfStrides=ep2.Stride_No .* sign(ep2.EarlyOrLate-.5);
    exemptLast=ep2.ExemptLast;
    exemptFirst=ep2.ExemptFirst;
    data=nan(length(labels),length(ep2));
    validStrides=nan(length(ep2),1);
    summaryFlag=ep2.summaryMethod;
    allData=cell(length(normalizedGroupData.ID),1);
    for s=1: length(normalizedGroupData.ID)
        
        dataPoints=normalizedGroupData.adaptData{s}.getEarlyLateData_v2(labels,conds,removeBiasFlag,numberOfStrides,exemptLast,exemptFirst,padWithNaNFlag);%Get data
        
        aux1(:,:,s)=squeeze(dataPoints{1});
    end
    
    % How we do it in the plotcheckerboards
    dataAll{l}=aux1;
    aux2=nanmean(aux1,1); %average across steps
    aux2=nanmedian(aux2,3); %median across participants
    aux2=reshape(aux2,12,28);
    [flippedEMGData] =flipEMGdata(aux2,1,2); %temporal of the fast leg
    aux2=flippedEMGData-fftshift(flippedEMGData,2);
    Yasym{l}=aux2(:,1:size(aux2,2)/2);
    Regressors{l}=Yasym{l};
    
     % How we do it in model
    dataAll{l}=aux1;
    aux=reshape(dataAll{l},size(dataAll{l},1),12,28,size(dataAll{l},3));
    aux=reshape(flipEMGdata(aux,2,3),size(aux,1),numel(labels),size(aux,4));
    aux2=nanmedian(aux,3); %median across participants
    aux2=aux2-fftshift(aux2,2);
    aux2=aux2(:,1:size(aux2,2)/2);
    aux3=nanmean(aux2,1); %average across steps
    Regressors2{l}=reshape(aux3,12,14);

    
    
    
    
end
%%
% figure 
% ytl={'GLU','TFL','HIP','RF','VL','VM','BF', 'SEMB','SEMT','MG','LG','SOL','PER','TA'};
% 
% % ytl=newLabelPrefix;
% % ytl=ytl(end:-1:1);
% yt=1:14;
% fs=14;
% 
% subplot(1,4,1)
% imagesc(Regressors{1}(:,end:-1:1)')
% caxis([-1 1])
% set(gca,'XTick',[],'YTick',yt,'YTickLabel',ytl,'FontSize',fs)
% title('Adaptation_{early}')
% 
% 
% 
% subplot(1,4,2)
% imagesc(Regressors{2}(:,end:-1:1)')
% caxis([-1 1])
% set(gca,'XTick',[],'YTick',yt,'YTickLabel',ytl,'FontSize',fs)
% title('Adaptation_{late}')
% colormap(flipud(map))
% 
% figure 
% ytl={'GLU','TFL','HIP','RF','VL','VM','BF', 'SEMB','SEMT','MG','LG','SOL','PER','TA'};
% yt=1:14;
% fs=14;
% 
% subplot(1,4,1)
% imagesc(Regressors2{1}(:,end:-1:1)')
% caxis([-1 1])
% set(gca,'XTick',[],'YTick',yt,'YTickLabel',ytl,'FontSize',fs)
% title('Adaptation_{early}')
% % title('F(A)')
% 
% 
% subplot(1,4,2)
% imagesc(Regressors2{2}(:,end:-1:1)')
% caxis([-1 1])
% set(gca,'XTick',[],'YTick',yt,'YTickLabel',ytl,'FontSize',fs)
% title('Adaptation_{late}')
colormap(flipud(map))
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
colormap(flipud(map))