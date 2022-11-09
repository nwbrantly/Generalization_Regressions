% This is a script to get the confidance interval of the step-by-step
% weights 

%% Load data and Plot checkerboard for all conditions.
% clear all; close all; clc;
clear all; clc;


% set script parameters, SHOULD CHANGE/CHECK THIS EVERY TIME.
groupID = 'BATS';
files = dir ([ groupID '*params.mat']);


n_subjects = size(files,1);

ii=0;
for i =1:n_subjects
    ii=1+ii;
    sub{ii} = files(i).name;
    subID{ii} = sub{ii}(1:end-10);
end

subID
n_subjects = size(files,1);

ep=defineRegressorsDynamicsFeedback('nanmean');
refEpTM = defineReferenceEpoch('TM base',ep);


%% Data Normalization 

mOrder={'TA','PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT', 'VM', 'VL', 'RF', 'HIP','TFL', 'GLU'};
GroupData=adaptationData.createGroupAdaptData(sub); %loading the data
GroupData=GroupData.removeBadStrides; %Removing bad strides
newLabelPrefix = defineMuscleList(mOrder); 
normalizedGroupData = GroupData.normalizeToBaselineEpoch(newLabelPrefix,refEpTM); %Normalized by OG base same as nimbus data
ll=normalizedGroupData.adaptData{1}.data.getLabelsThatMatch('^Norm');
l2=regexprep(regexprep(ll,'^Norm',''),'_s','s');
normalizedGroupData=normalizedGroupData.renameParams(ll,l2);
newLabelPrefix = regexprep(newLabelPrefix,'_s','s');



%%  Getting the step-by-step data 

%Adaptation epochs
strides=[-40 440 200];
cond={'TM base','Adaptation','Post 1'}; %Conditions for this group


exemptFirst=[1];
exemptLast=[5]; %Strides needed
names={};
shortNames={};

ep=defineEpochs(cond,cond,strides,exemptFirst,exemptLast,'nanmean',{'Base','Adapt','Post1'}); %epochs

padWithNaNFlag=true;
[dataEMG,labels,allDataEMG]=normalizedGroupData.getPrefixedEpochData(newLabelPrefix(end:-1:1),ep,padWithNaNFlag);
%Flipping EMG:
for i=1:length(allDataEMG)
    aux=reshape(allDataEMG{i},size(allDataEMG{i},1),size(labels,1),size(labels,2),size(allDataEMG{i},3));
    allDataEMG{i}=reshape(flipEMGdata(aux,2,3),size(aux,1),numel(labels),size(aux,4));
end

EMGdata=cell2mat(allDataEMG);
muscPhaseIdx=1:size(EMGdata,2);%336; %14 muscles per leg %360; %All muscles
%% Getting the C values 
epochOfInterest={'TM base','TM mid 1','PosShort_{early}','PosShort_{late}','Ramp','Optimal','Adaptation','Adaptation_{early}','TiedPostPos','TMmid2','NegShort_{late}','Post1_{Early}','TMbase_{early}'};
ep=defineRegressorsDynamicsFeedback('nanmean');
refEpTM = defineReferenceEpoch('TM base',ep);
flip=1;

if flip==1
n=2;
method='IndvLegs';
else
   n=1; 
   method='Asym';
end

for s=1:n_subjects
    for l=1:length(epochOfInterest)
        ep2=defineReferenceEpoch(epochOfInterest{l},ep);
        adaptDataSubject = normalizedGroupData.adaptData{1, s};
        [~,~,~,Data{s,l}]=adaptDataSubject.getCheckerboardsData(newLabelPrefix,ep2,[],flip);
   end
end

%% Group C

figure(1)
hold on
bootstrap=1
X1=[];
X2=[];
replacement=0





temp=[];
temp2=[];
DataBoot={};


DataBoot=Data(1:n_subjects,:);

for c=1:length(Data)
    for s=1:n_subjects
        temp(:,:,s)=DataBoot{s,c};
    end
    temp2{1,c}=nanmedian(temp,3);
    tt=temp2{1,c}(:,end:-1:1);
    temp2{1,c}=tt;
    temp2{1,c}=reshape(temp2{1,c},14*2*12,1);
end

temp2=cell2mat(temp2)';
temp2=temp2-fftshift(temp2,2);
temp2=temp2(:,1:size(temp2,2)/2,:);
temp2=temp2';
C2=[temp2(:,5) temp2(:,6)];

for i=1:n_subjects
    Y2=EMGdata(:,muscPhaseIdx,i);
    Y2=nanmedian(Y2,3);
    Yasym2=Y2-fftshift(Y2,2);
    Yasym2=Yasym2(:,1:size(Yasym2,2)/2,:);
    
    
    Cinv2=pinv(C2);
    Xdynamics= Cinv2*Yasym2'; %x= y/C
    
    X1=[X1;Xdynamics(1,:)];
    X2=[X2;Xdynamics(2,:)];
    
    subplot(2,1,1)
    hold on
    plot(Xdynamics(1,:))
    subplot(2,1,2)
    hold on
    plot(Xdynamics(2,:))
    
    
    
    
end


legend(subID,'AutoUpdate','off') 
    subplot(2,1,1)
    hold on
    yline(0)
    ylabel('W_{reactive}')
        subplot(2,1,2)
    hold on
    yline(0)
    ylabel('W_{context}')

set(gcf,'color','w')

%% Individual fit with individual C


figure(1)
hold on
bootstrap=1
X1=[];
X2=[];
replacement=0






temp=[];

DataBoot={};

for i=1:n_subjects
    temp2=[];
    
    DataBoot=Data(1:n_subjects,:);
    
    for c=1:length(Data)
        temp2{1,c}=DataBoot{i,c};
        tt=temp2{1,c}(:,end:-1:1);
        temp2{1,c}=tt;
        temp2{1,c}=reshape(temp2{1,c},14*2*12,1);
    end
    
    temp2=cell2mat(temp2)';
    temp2=temp2-fftshift(temp2,2);
    temp2=temp2(:,1:size(temp2,2)/2,:);
    temp2=temp2';
    C2=[temp2(:,5) temp2(:,6)];
    
    
    Y2=EMGdata(:,muscPhaseIdx,i);
    Y2=nanmedian(Y2,3);
    Yasym2=Y2-fftshift(Y2,2);
    Yasym2=Yasym2(:,1:size(Yasym2,2)/2,:);
    
    
    Cinv2=pinv(C2);
    Xdynamics= Cinv2*Yasym2'; %x= y/C
    
    X1=[X1;Xdynamics(1,:)];
    X2=[X2;Xdynamics(2,:)];
    
    subplot(2,1,1)
    hold on
    plot(Xdynamics(1,:))
    subplot(2,1,2)
    hold on
    plot(Xdynamics(2,:))
    
    
    
    
end

legend(subID,'AutoUpdate','off') 
    subplot(2,1,1)
    hold on
    yline(0)
    ylabel('W_{reactive}')
        subplot(2,1,2)
    hold on
    yline(0)
    ylabel('W_{context}')

set(gcf,'color','w')


%%

grayColor = [.7 .7 .7];

figure(2)
subplot(2,1,1)
hold on
% X_1=[X1(:,1:40) nan(size(X1(:,1:40),1),1) X1(:,41:480) nan(size(X1(:,1:40),1),1) X1(:,481:end)];
% y = nanmean(X1,1); % your mean vector;
% x = 1:numel(y);
% std_dev = nanstd(X1,1);
% curve1 = y + std_dev;
% curve2 = y - std_dev;
% x2 = [x, fliplr(x)];
% inBetween = [curve1, fliplr(curve2)];
% fill(x2, inBetween,'r','FaceAlpha',0.3,'EdgeColor','none')
% hold on;
% ylabel('W_{reactive}')
% xlabel('Strides')
% plot(x, y, 'r', 'LineWidth', 2);

d{1} = X1(:,1:40);
d{2} = X1(:,41:480);
d{3} =  X1(:,481:end);
index{1}= 1:40;
index{2}= 42:481;
index{3}= 483:682;

for i=1:3
    y = nanmean(d{i},1); % your mean vector;
    %     x = 1:numel(y);
    x=index{i};
    std_dev = nanstd(d{i},1);
    curve1 = y + std_dev;
    curve2 = y - std_dev;
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    fill(x2, inBetween,'r','FaceAlpha',0.3,'EdgeColor','none')
    hold on;
    ylabel('W_{reactive}')
    xlabel('Strides')
    plot(x, y, 'r', 'LineWidth', 2);
   
end
yline(0)
% legend('Training')


subplot(2,1,2)
hold on
% X_2=[X2(:,1:40) nan(size(X2(:,1:40),1),1) X2(:,41:480) nan(size(X2(:,1:40),1),1) X2(:,481:end)];
% y = nanmean(X2,1); % your mean vector;
% x = 1:numel(y);
% std_dev = nanstd(X2,1);
% curve1 = y + std_dev;
% curve2 = y - std_dev;
% x2 = [x, fliplr(x)];
% inBetween = [curve1, fliplr(curve2)];
% fill(x2, inBetween,'r','FaceAlpha',.3,'EdgeColor','none')
% hold on;
% ylabel('W_{context}')
% xlabel('Strides')
% plot(x, y, 'r', 'LineWidth', 2);
d{1} = X2(:,1:40);
d{2} = X2(:,41:480);
d{3} =  X2(:,481:end);

for i=1:3
    y = nanmean(d{i},1); % your mean vector;
    %     x = 1:numel(y);
    x=index{i};
    std_dev = nanstd(d{i},1);
    curve1 = y + std_dev;
    curve2 = y - std_dev;
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    fill(x2, inBetween,'r','FaceAlpha',0.3,'EdgeColor','none')
    hold on;
    ylabel('W_{context}')
    xlabel('Strides')
    plot(x, y, 'r', 'LineWidth', 2);
    
end
yline(0)
set(gcf,'color','w')
% save([groupID,'_',num2str(n_subjects),'_iteration_', num2str(n)],'X1','X2')
%%

disp('BATS - Reactive')
[h,p,~,~] = ttest(X1(:,482))
disp('BATS - Context')
[h,p,~,~] = ttest(X2(:,482))
xlabel('Context')

%%


disp('BATR - Reactive')
[h,p,~,~] = ttest(X1(:,482))
disp('BATS - Context')
[h,p,~,~] = ttest(X2(:,482))


% set(gcf,'color','w')