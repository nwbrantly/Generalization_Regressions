%Model fit per muscle

%% Adding path 

% upload your path
main='/Users/dulcemariscal/Documents/GitHub/';
% main='C:\Users\dum5\OneDrive - University of Pittsburgh\_BoxMigration\GitHub\';
addpath(genpath([main, 'Generalization_Regressions']))
addpath(genpath([main,'labTools']))
% addpath(genpath([main,'LongAdaptation']))
addpath(genpath([main,'splitbelt-EMG-adaptation']))
addpath(genpath([main,'EMG-LTI-SSM']))
addpath(genpath([main,'matlab-linsys']))
% rmpath(genpath([main,'PittSMLlab']))

%% Load real data:
clear all;clc;close all
%% Free model - Linear regression - Asymmetry with baseline
%% This is just the saved data - Update accrodingly 


groupID='C3'

%EMG Data 
% fname = ['dynamicsData_', groupID, '_subj_3_RemoveBadMuscles_0.h5']
fname='dynamicsData_C3_subj_3_RemoveBadMuscles_0.h5'

%%%%%%%% Regressors
% load([ groupID,'_1_IndvLegsC4_RemovedBadMuscle_0RemoveBias_1.mat'])


%%
muscleOrder={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT',...
    'VM', 'VL', 'RF', 'HIP','TFL', 'GLU'};
ytl=([strcat('f',muscleOrder) strcat('s',muscleOrder)]);  %List of muscle
ytl(end:-1:1) = ytl(:);
yt=1:length(ytl);
fs=14;
%% Getting the data

EMGdata=h5read(fname,'/EMGdata');
binwith=10;
[Y,~,U,Ubreaks,Ysum,Yindv,labels,C,Cinv]=groupDataToMatrixForm_Update(1:size(EMGdata,3),fname,0);

% color=colormap(jet(size(Yindv,3)));
poster_colors;
colorOrder=[p_red; p_orange; p_fade_green; p_fade_blue; p_plum; p_green; p_blue; p_fade_red; p_lime; p_yellow; [0 0 0];[0 1 1];p_red; p_orange; p_fade_green; p_fade_blue; p_plum; p_green; p_blue; p_fade_red; p_lime; p_yellow; [0 0 0];[0 1 1]];
c=[ colorOrder; colorOrder;colorOrder];
subjecname={'C3S01','C3S02','C3S03'};


%Pick the conditions that we are going to use for the regression
epochOfInterest=h5read(fname,'/Epochs')';
epochOfInterest=deblank(epochOfInterest); %Remove random null values added during the saving of the matrix

context= find(strcmp(strtrim(epochOfInterest) ,'Adaptation')==1);
reactive=find(strcmp(epochOfInterest,'NegShort_{early}')==1);
tmbase=find(strcmp(epochOfInterest,'TM base')==1);

for s=1:size(Yindv,3)
Yasym=Yindv(:,:,s);
C=Cinv(:,:,s)';
Uf=[U;ones(size(U))];

%% Organizing the data

removebaseline=1; %Remove bias flag


Casym=[C(:,reactive) C(:,context)]; % EMGreactive and EMGcontext
Ymodel=Yasym';  % Transpose the EMG data to match equations


if removebaseline==1
    bias=nanmean(Yasym(1:35,:)); %Computing the bias
    Casym=Casym-C(:,tmbase);
%     Casym=Casym-bias'; %Removing bias from the regressors
    Ymodel=Ymodel-bias'; %removing bias from the data
end

    %Organizing the data per muscle
    Cmuscles=reshape(Casym',2,12,size(Yasym,2)/12); %muscles for regressors 2 number of regressors, 12 number of phase of the gait cycle
    Ymuscles=reshape(Ymodel(:,1:size(Yasym,1))',size(Yasym,1),12,size(Yasym,2)/12); %data



%% Linear regression individual muscles
reconstruction_indv=[];
data=[];
C_indv=[];
X_indv=[];

for i=1:28
    unit(:,:,i)=Cmuscles(:,:,i)'./vecnorm(Cmuscles(:,:,i)'); %Getting the unit vector of the regressors
    %         temp(:,:,i)=pinv(unit(:,:,i)'); % Getting the inverse
        if sum(any(isnan(unit(:,:,i))))>=1
            temp(:,:,i)=(unit(:,:,i)); %if nan can't get inverse, carry the nan alone.
        else
            temp(:,:,i)=pinv(unit(:,:,i)'); % Getting the inverse, if data is not nan.
        end
        
    X2asym(:,:,i) =temp(:,:,i)'*Ymuscles(:,:,i)'; %x= y/C
    Y2asym(:,:,i)=  unit(:,:,i)* X2asym(:,:,i) ; %yhat
    temp2(:,:,i)=X2asym(:,:,i)'; % transposing the dynamics vector
    model{i}.C=unit(:,:,i); %saving the regressors. Yhis is necesary for the plotting funciton
    
    reconstruction_indv =[ reconstruction_indv ; Y2asym(:,:,i)]; % Concatenating the data reconstructed
    data =[ data ; Ymuscles(:,:,i)'];  % Concatenating the data
    C_indv=[C_indv;unit(:,:,i)];  % Concatenating the regressors for each muscles
    X_indv=[X_indv,X2asym(:,:,i) ];  % Concatenating the dynamics
    
    % Checking for colinearity and correlation between the regresssions
%     temp5=corrcoef(model{i}.C); % correlation coeficient
%     correlation(i,1)=temp5(2);
%     temp10(i,:)=vif([model{i}.C nanmean(Ymuscles(481:491,:,i))']); %Variance inflation
%     temp5=vif([model{i}.C]);
%     impact(i,:)=temp5;
    
    %fit linear model 
    mdl{i}= fitlm(unit(:,:,i),nanmean(Ymuscles(41:45,:,i))','VarNames',{'Reactive','Contextual', labels(i).Data(1:end-1)},'Intercept',false);
    mdl2{i}= fitlm(unit(:,:,i),nanmean(Ymuscles(41:45,:,i))','VarNames',{'Reactive','Contextual', labels(i).Data(1:end-1)},'Intercept',true);
   
    figure(1)
    hold on
    li{s}=scatter(i,squeeze(nanmean(temp2(41:45,1,i))),100,"filled",'MarkerFaceColor', c(s,:));
    errorbar(i,mdl{i}.Coefficients.Estimate(1),mdl{i}.Coefficients.SE(1),'k')
    
    figure(2)
    hold on
    li2{s}=scatter(i,squeeze(nanmean(temp2(41:45,2,i))),100,"filled",'MarkerFaceColor', c(s,:));
    errorbar(i,mdl{i}.Coefficients.Estimate(2),mdl{i}.Coefficients.SE(2),'k')
    
    
end

   
    data=nanmean(data(:,41:45),2);
    data = data(~isnan(data));
    
    reconstruction_indv=nanmean(reconstruction_indv(:,41:45),2);
    reconstruction_indv=reconstruction_indv(~isnan(reconstruction_indv));
    
    reconstruction_reactive=nanmean(reconstruction_reactive(:,41:45),2);
    reconstruction_reactive=reconstruction_reactive(~isnan(reconstruction_reactive));
    
    
    r2_total(subj,1)=my_Rsquared_coeff(data,reconstruction_indv);
    r2_reactive(subj,1)=my_Rsquared_coeff(data,reconstruction_reactive);
    
    r2_total_uncenter(subj,1)=Rsquared_uncenter(data,reconstruction_indv);
    r2_reactive_uncenter(subj,1)=Rsquared_uncenter(data,reconstruction_reactive);
    

% index=find(impact>5); % fidn the regressor with high colinearity

Uf=Uf(:,1:size(temp2,1));
set(gcf,'color','w')
Uf=Uf(:,1:size(temp2,1));

%     
%     if size(Yindv,3)>1
%     
%         figure(1)
%         hold on
%         scatter(1:28,squeeze(nanmean(temp2(41:45,1,:))),100,"filled");
%         
%         
%         figure(2)
%         hold on 
%         scatter(1:28,squeeze(nanmean(temp2(41:45,2,:))),100,"filled");
%         
        if s==size(Yindv,3)
            figure(1)
            ylabel({'Reactive';'AU'})
            yline(0)
            set(gca,'XTick',yt,'XTickLabel',ytl,'FontSize',10)
            legend([li{:}],subjecname, 'location','Best')
            set(gcf,'color','w')
            
            figure(2)
            ylabel({'Contextual';'AU'})
            yline(0)
            set(gca,'XTick',yt,'XTickLabel',ytl,'FontSize',10)
            legend([li2{:}],{'C3S01','C3S02','C3S03'}, 'location','Best')
            set(gcf,'color','w')
        end
        
    
    
%     end
    


end
%%  Plotting 
model2{1}.C=C_indv; %Save regressors for model
model2{1}.Out= reconstruction_indv; %Save reconstruction of the indiivudal muscles to the model format for plotting 
analysis=0; %Flag asking if you want to run the analysis
isF=0; %Flag for fast leg


legacy_vizSingleModel_FreeModel_ShortAdaptation_IndvLeg2(model2{1},data,Uf,analysis,[],isF)
% legacy_vizSingleModel_FreeModel_ShortAdaptation_IndvLeg(model2{1},data,Uf,analysis,[],isF)
% save('BATS_indv_muscles.mat','X2asym')


%% Plotting time course for individual muscles
analysis=0

for i=1:28
%pick the muscle that you want
    
    % Pick the data that you want to plot 
    %    Xasym=[temp2(1:40,:,i);nan(1,size(temp2,2));temp2(41:480,:,i);nan(1,size(temp2,2));temp2(481:end,:,i)];
    Xasym=[temp2(41:490,:,i)];
    %      Xasym=[temp2(1:200,:,i)];
    % Xasym=[X2asym(1:40,:);nan(1,size(X2asym,2));X2asym(41:80,:);...
    %     nan(1,size(X2asym,2));X2asym(81:520,:);nan(1,size(X2asym,2));X2asym(521:end,:)];
    % Xasym=[X2asym(681:843,:);nan(1,size(X2asym,2))];
    
    figure
    subplot(2,1,1)
    hold on
    % scatter(1:length(movmean(Xasym(:,1),binwith)), movmean(Xasym(:,1),binwith),10,'k','filled')
    scatter(1:length(movmean(Xasym(:,1),binwith)), movmean(Xasym(:,1),binwith),'filled','MarkerFaceColor',"#EDB120") %"#77AC30" )%
        
    % legend('Baseline','AutoUpdate','off')
    legend('Negative')
    title(labels(i).Data)
    %     legend('Removal Perturbation','AutoUpdate','off')
    % uistack(pp,'bottom')
    yline(0)
    ylabel({'Reactive';'(A.U)'})
    %     ylabel({'Removal';'(A.U)'})
    xlabel('strides')
    
    if size(temp2(:,:,i),2)>=2
        % figure
        subplot(2,1,2)
        hold on
        scatter(1:length(movmean(Xasym(:,2),binwith)), movmean(Xasym(:,2),binwith),'filled','MarkerFaceColor'," #00008B")
        %
        legend('Contextual')
        % % legend('Switch','AutoUpdate','off')
        % % uistack(pp,'bottom')
        yline(0)
        ylabel({'Contextual';'(A.U)'})
        xlabel('strides')
    end
    if i<=14
        isF=0;
    else
        isF=1;
    end
    set(gcf,'color','w')
    % Plot the time course plus the regressors 
%     legacy_vizSingleModel_FreeModel_ShortAdaptation_IndvLeg(model{i},Ymuscles(:,:,i)',Uf,analysis,{labels(i).Data(2:end-1)},isF)
    set(gcf,'color','w')

end



%% 
%% Geting the average of the first 5 strides post-adaptation per muscle
figure
hold on
for i=1:28
    if i<15
        Li{1}=scatter(nanmean(temp2(41:45,1,i)),nanmean(temp2(41:45,2,i)),100,"filled",'MarkerFaceColor', 'b');
        text(nanmean(temp2(41:45,1,i)),nanmean(temp2(41:45,2,i)),{labels(i).Data(2:end-1)})
    else
        Li{2}=scatter(nanmean(temp2(41:45,1,i)),nanmean(temp2(41:45,2,i)),100,"filled",'MarkerFaceColor', 'r')  ;
        text(nanmean(temp2(41:45,1,i)),nanmean(temp2(41:45,2,i)),{labels(i).Data(2:end-1)})
    end
    
end
legend([Li{:}],['Slow';'Fast'])
ylabel({'Contextual';'A.U'})
xlabel({'Reactive';'A.U'})

% xlim([-0.3 1.4])
% ylim([-0.3 1.4])
set(gcf,'color','w')

%% Geeting the average of Early Adapt and lateAdapt
figure
hold on
for i=1:28
    if i<15
        Li{1}=scatter(nanmean(temp2(41:45,1,i)),nanmean(temp2(41:45,2,i)),100,"filled",'MarkerFaceColor', 'b');
        text(nanmean(temp2(41:45,1,i)),nanmean(temp2(41:45,2,i)),{labels(i).Data(2:end-1)})
    else
        Li{2}=scatter(nanmean(temp2(41:45,1,i)),nanmean(temp2(41:45,2,i)),100,"filled",'MarkerFaceColor', 'r')  ;
        text(nanmean(temp2(41:45,1,i)),nanmean(temp2(41:45,2,i)),{labels(i).Data(2:end-1)})
    end
    
end
legend([Li{:}],['Slow';'Fast'])
ylabel({'Contextual';'A.U'})
xlabel({'Reactive';'A.U'})
% axis square
xlim([-1 3])
ylim([-1 1])

title('EarlyAdapt')
set(gcf,'color','w')


figure
hold on
for i=1:28
    if i<15
        Li{1}=scatter(nanmean(temp2(440:480,1,i)),nanmean(temp2(440:480,2,i)),100,"filled",'MarkerFaceColor', 'b');
        text(nanmean(temp2(440:480,1,i)),nanmean(temp2(440:480,2,i)),{labels(i).Data(2:end-1)})
    else
        Li{2}=scatter(nanmean(temp2(440:480,1,i)),nanmean(temp2(440:480,2,i)),100,"filled",'MarkerFaceColor', 'r')  ;
        text(nanmean(temp2(440:480,1,i)),nanmean(temp2(440:480,2,i)),{labels(i).Data(2:end-1)})
    end
    
end
legend([Li{:}],['Slow';'Fast'])
ylabel({'Contextual';'A.U'})
xlabel({'Reactive';'A.U'})
axis square
title('LateAdapt')
set(gcf,'color','w')
xlim([-.1 .85])
ylim([-.1 .85])

%% Comparing individual analysis with the leg specific analysis

Cslow=[Casym(1:size(C,1)/2,1)  Casym(1:size(C,1)/2,2)]./vecnorm([Casym(1:size(C,1)/2,1)  Casym(1:size(C,1)/2,2)]); %SLOW
Cfast=[Casym(1+size(C,1)/2:end,1)  Casym(1+size(C,1)/2:end,2)]./vecnorm([Casym(1+size(C,1)/2:end,1)  Casym(1+size(C,1)/2:end,2)]); %FAST

Ymodel=Ymodel';
Yslow=Ymodel(1:680,1:size(Ymodel,2)/2); %SLOW
Yfast=Ymodel(1:680,size(Ymodel,2)/2+1:end); %FAST



Cinv=pinv(Cslow);
Xslow = Cinv*Yslow'; %x= y/C
hatYslow=  Cslow *Xslow ; %yhat = C

% hatYslow= hatYslow';

Cinv=pinv(Cfast);
Xfast = Cinv*Yfast'; %x= y/C
hatYfast=  Cfast *Xfast ; %yhat = C

% hatYfast=hatYfast;
Yfast=Yfast';
Yslow=Yslow';

% Variance explained

ex2=[0.2314    0.2980    0.7529];
ex1=[0.7255    0.0863    0.1608];
mid=ones(1,3);
N=100;
gamma=1.5; %gamma > 1 expands the white (mid) part of the map, 'hiding' low values. Gamma<1 does the opposite
gamma=1;
map=[flipud(mid+ (ex1-mid).*([1:N]'/N).^gamma); mid; (mid+ (ex2-mid).*([1:N]'/N).^gamma)];



%% Variance all muscles
Rsquared_slow2= my_Rsquared_coeff(Yslow,hatYslow);
Rsquared_fast2 = my_Rsquared_coeff(Yfast,hatYfast);

Rsquared_slow_unit2 = my_Rsquared_coeff(Yslow,reconstruction_indv(1:168,:));
Rsquared_fast_unit2 = my_Rsquared_coeff(Yfast,reconstruction_indv(169:end,:));

figure()
plot(Rsquared_slow2)
hold on
plot(Rsquared_slow_unit2)
legend('Per leg','Per muscle')

figure
plot(Rsquared_fast2)
hold on
plot(Rsquared_fast_unit2)
legend('Per leg','Per muscle')

% both legs
Rsquared_both= my_Rsquared_coeff([Yslow(:,:);Yfast(:,:)],[hatYslow(:,:);hatYfast(:,:)]);
Rsquared_both_unit = my_Rsquared_coeff([Yslow(:,:);Yfast(:,:)],reconstruction_indv(:,:));

figure
plot(Rsquared_both)
hold on
plot(Rsquared_both_unit)
legend('Per leg','Per muscle')
title('Both leg')
set(gcf,'color','w')
% figure(3)
% ytl={'GLU','TFL','HIP','RF','VL','VM','BF', 'SEMB','SEMT','MG','LG','SOL','PER','TA'};
%%
%% Calculating the R squared per muscle
% mm= 0:12:168;
% mm2=1:12:168;
% for m=1:14
%     Rsquared_slow(m,:) = my_Rsquared_coeff(Yslow(mm2(m):mm(m+1),:),hatYslow(mm2(m):mm(m+1),:));
%     Rsquared_fast(m,:) = my_Rsquared_coeff(Yfast(mm2(m):mm(m+1),:),hatYfast(mm2(m):mm(m+1),:));
%     
%     Rsquared_slow_unit(m,:) = my_Rsquared_coeff(Ymuscles(:,:,m)',Y2asym(:,:,m));
%     Rsquared_fast_unit(m,:) = my_Rsquared_coeff(Ymuscles(:,:,m+14)',Y2asym(:,:,m+14));
%     
% end
% % figure
% muscleOrder={'GLU','TFL','HIP','RF','VL','VM','BF', 'SEMB','SEMT','MG','LG','SOL','PER','TA'};
% ytl= defineMuscleListV2(muscleOrder); %List of muscle
% % ytl={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT', 'VM', 'VL', 'RF', 'HIP','TFL', 'GLU'};
% % ytl(end:-1:1) = ytl(:);
% binw=5;
% 
% mtp=10;
% for m=mtp
%     %     figure
%     %      imagesc(nanmean(Yslow(mm2(m):mm(m+1),481:485),2)',[-1 ,1])
%     %      colormap(flipud(map))
%     figure
%     %     subplot(14,1,m)
%     hold on
%     aux1=conv(Rsquared_slow(m,:),ones(1,binw)/binw,'valid'); %Smoothing
%     plot(aux1,'LineWidth',2,'DisplayName','Per leg','Color',"#0072BD") ;
%     
%     aux1=conv(Rsquared_slow_unit(m,:),ones(1,binw)/binw,'valid'); %Smoothing
%     plot(aux1,'LineWidth',2,'DisplayName','Per muscle','Color',"#A2142F") ;
%     
%     %     plot(movmean(muscles_slow(m,:),5))
%     %     plot(movmean(muscles_slow_unit(m,:),5))
%     %     scatter(1:length(muscles(m,:)),movmean(muscles(m,:),5),'filled')
%     %     legend({'Per leg';'Per muscle'},'AutoUpdate','off');
%     legend('Location','NorthEastOutside','AutoUpdate','off')
%     ylabel(ytl{m+14})
%     %     yline(nanmean(muscles_slow(m,10:30)))
%     
%     pp=patch([40 480 480 40],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
%     
%     %     pp=patch([0 440  440 0],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
%     
%     uistack(pp,'bottom')
% end
% set(gcf,'color','w')
% 
% % figure
% for m=mtp
%     %     subplot(14,1,m)
%     figure
%     hold on
%     aux1=conv(Rsquared_fast(m,:),ones(1,binw)/binw,'valid'); %Smoothing
%     plot(aux1,'LineWidth',2,'DisplayName','Per leg','Color',"#0072BD") ;
%     
%     aux1=conv(Rsquared_fast_unit(m,:),ones(1,binw)/binw,'valid'); %Smoothing
%     plot(aux1,'LineWidth',2,'DisplayName','Per muscle','Color',"#A2142F") ;
%     %         plot(movmean(muscles_fast(m,:),5))
%     %         plot(movmean(muscles_fast_unit(m,:),5))
%     %     scatter(1:length(muscles(m,:)),movmean(muscles(m,:),5),'filled')
%     %     legend({'Per leg';'Per muscle'},'AutoUpdate','off');
%     
%     legend('Location','NorthEastOutside','AutoUpdate','off')
%     ylabel(ytl{m})
%     %     yline(nanmean(muscles_fast(m,10:30)))
%     pp=patch([40 480 480 40],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
%     %     pp=patch([0 440  440 0],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
%     
%     uistack(pp,'bottom')
% end
% set(gcf,'color','w')