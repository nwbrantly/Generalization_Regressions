%% Script to plot the reactive and contexutal recruitment time courses and the model R2 
%% This script needs the input that it is save from model_fit_individual_muscle_PerSubject.m 

% Define which group you want to plot
% groupID={'VATR','VATS'};
% groupID={'NTR','NTS'};
%  groupID={'CTR','CTS'};
groupID={'C3S10'}; 

f=[1:14 1:14;1:14 1:14]; %Hard coded to make it look cute! 

binwith=5; % Running windown length 

for g=1:length(groupID) %Group loop (Device or W/O device)
    
    % Load data
    if strncmpi(groupID{g},'BAT',3) %See if the first letters of the group are BAT
        load([groupID{g},'_post-adaptation_Indv_0_04-April-2024'])
    else
        load([groupID{g}, '_post-adaptation_Indv_0_22-April-2024.mat'])
    end
        
    for muscle=1:length(labels) %% Define muscles to plot (This data set has a total of 28 muscles)
        
        r2 = my_Rsquared_coeff(model{muscle}.EMGobserved', model{muscle}.EMGestimated',1); % Compute the R2 throuhgtout the data
        figure(f(g,muscle))
        
        % Pick the data that you want to plot
        Wdata(:,1)=[reactive_trace(41:end,muscle)]'; %Data organization and selecting only the postadaptation data
        Wdata(:,2)=[contextual_trace(41:end,muscle)]'; %Data organization and selecting only the postadaptation data
        r2=[r2(:,41:end)];
        
        % Making indeces where r2 is negative NaN
        idx= find(r2<0); % finding when the r2 is negative
        Wdata(idx,1)=nan; %reactive
        Wdata(idx,2)=nan; %contextual
        r2(:,idx)=nan; %r2 
        
        % Getting the average of the data while ignoting nan. Note that
        % overrepresentation of data can happen
        avg(:,1)=movmean(Wdata(:,1),binwith,'omitnan'); %reactive
        avg(:,2)=movmean(Wdata(:,2),binwith,'omitnan'); %contextual
        r2_avg= movmean(r2,binwith,'omitnan');
        
        % Alternative: we get the average and then we mark data as NaN.
        % This leads to more empty data points
        %             idx= find(r2_avg<0);
        %             avg(idx,:)=nan;
        %             r2_avg(:,idx)=nan;
        
        
        %Plot reactive time course 
        if muscle<15 %Choosing where to plot the time course 
            subplot(size(Wdata,2)+1,2,1)
        else
            subplot(size(Wdata,2)+1,2,2)
        end
        
        hold on
        scatter(1:length(avg),avg(:,1),'filled','DisplayName',groupID{g}); %plot  
        legend %Adding the groupID 
        title(labels(muscle).Data) %adding muscle that we are looking at on the title 
        yline(0,'HandleVisibility','off')
        ylabel({'Reactive';'(A.U)'})
        xlabel('strides')
              
        %Plot  contextual time course 
        if size(Wdata,2)>=2
            if muscle<15
                subplot(size(Wdata,2)+1,2,3)
            else
                subplot(size(Wdata,2)+1,2,4)
            end
            hold on
            scatter(1:length(avg),avg(:,2),'filled','DisplayName',groupID{g}); % plot 
            legend  %Adding the groupID 
            yline(0,'HandleVisibility','off')
            ylabel({'Contextual';'(A.U)'})
            xlabel('strides')
        end
        
        % Plot the time course of the R2
        if muscle<15
            subplot(size(Wdata,2)+1,2,5)
        else
            subplot(size(Wdata,2)+1,2,6)
        end
        
        hold on
        plot(1:length(r2_avg),r2_avg,'DisplayName',groupID{g}); %plot 
        ylim([-.5 1])
        yline(0,'HandleVisibility','off')
        legend
        xlabel('strides')
        ylabel('R^2')
        set(gcf,'color','w')        
       
    end
   
end

%% Plot of the data and reconstruction per muscle 

groupID={'C3S10'}; % Define which group you want to plot
%Define type analysis: 
%0 or else: Linear regression with input regressors 
%1:Reactive and adaptation without dynamics 
%2:Analysis using PCA in all the data 
%3: removing the mean on the data. similar analysis as PCA     

analysis=0; 

% Load data
if strncmpi(groupID{g},'BAT',3) %See if the first letters of the group are BAT
    load([groupID{g},'_post-adaptation_Indv_0_04-April-2024'])
else
    load([groupID{g},'_post-adaptation_Indv_0_22-April-2024.mat'])
end

for muscle=1  %% Define muscles to plot (This data set has a total of 28 muscles)
    
    if muscle<=14
        isF=0; % Slow leg
    else 
        isF=1; % Fast leg
    end
    
    % Plot the time course plus the regressors
    vizSingleModel_LinearRegression_IndvLeg(model{muscle},model{muscle}.EMGobserved',analysis,{labels(muscle).Data(1:end-1)},isF)
    
end

