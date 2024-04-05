clear all
poster_colors;
colorOrder=[p_red; p_orange; p_fade_green; p_fade_blue; p_plum; p_green; p_blue; p_fade_red; p_lime; p_yellow; [0 0 0]];

sqrtFlag=false;
loop=0;
Li=[];

labels={'YA_{TR}','YA_{TS}','YA'};
groups=[2];
plotindv=1;
timeCourse=1;
binwith=5;

%  figure(1)
for group= groups
    loop=loop+1;
    C=[];
    modelRed=[];
    
    if group==1
        fname='dynamicsData_ATR_V4.h5';
        subj=[1:4];
        %         load('YA_TR_fixDandCV4_20220316T114557.mat')
        %         load A_15_AsymC3_EarlyLateAdaptation.mat
        load ATR_4_AsymC3_EarlyLateAdaptation
        load('ATR_fixDandCV4_20220316T114557.mat')
        color=[0 0.4470 0.7410];
        
    elseif group==2
        %         fname='dynamicsData_ATS_V6.h5';
        subj=[1:9 11];
        load ATS_Fix_C&D140422T182527
        %                 load('YA_TS_fixCandD1_280322T212228.mat')
        %         load A_15_AsymC3_EarlyLateAdaptation.mat
        load ATS_11_AsymC3_EarlyLateAdaptation.mat
        color=[0.8500 0.3250 0.0980];
        %         color=p_fade_blue;
        
    elseif group==3
        fname='dynamicsData_ATall.h5';
        subj=[1:6 8:13];
        %         load('A_15_Asym_EarlyLate40Adaptation.mat')
        load A_15_AsymC3_EarlyLateAdaptation.mat
        color=p_fade_red;
    end
    
    X1=[];
    X2=[];
    i=0;
    for s=1
        
        
        
        %         i=i+1;
        %         subjIdx=[subj]; %subjects tgat we are using 
        %         [Y,Yasym,Ycom,U,Ubreaks]=groupDataToMatrixForm(subjIdx,sqrtFlag,fname);
        %         %getting data 
        %         
        %         Uf=[U;ones(size(U))];
        %         datSet=dset(Uf,Yasym');
        
        
        
        Yasym=datSet.out';
        
        %% Free model - Linear regression - Asymmetry
        modelRed.C=C(:,1:3);
        if ~isempty(modelRed)
            C=modelRed.C;
        end
        
        Cinv=pinv(C)';
        X2asym = Yasym*Cinv; %x= y/C
        Y2asym= C * X2asym' ; %yhat = C
        
        for i=1:size(X2asym,2)
           
           
            X1(:,i)=X2asym(:,1);
            
            
            X2(:,i)=X2asym(:,2);
            
            X3(:,i)=X2asym(:,3);
        end
        
        if plotindv
            figure(2)
            subplot(3,1,1)
            hold on
            plot( movmean(X2asym(:,1),binwith))
            
            subplot(3,1,2)
            plot(movmean(X2asym(:,2),binwith))
            hold on
            
            subplot(3,1,3)
            plot(movmean(X2asym(:,3),binwith))
            hold on
            
        end
        
        
    end
    %%
    if timeCourse
        figure(1)
        cond={'TMbase','Adaptation','AdaptationEnd','Post1'};
        if group==1 || group==2 || group==3 
            %         strides=[50 950 1140];
            %         ini=[1 51 951];
            %         sz=[50 900 190];
            strides=[50 900 950 1140];
            ini=[1 51 901 951];
            sz=[50 850 50 190];
        else
            strides=[50 600 800];
            ini=[1 51 601];
            sz=[50 550 200];
            
            %         strides=[50 600 800];
            %         ini=[1 51 601];
            %         sz=[50 550 200];
            
            strides=[50 550 600 800];
            ini=[1 51  551 601];
            sz=[50 500 50 200];
            
            
        end
        Opacity=0.5;
        %
        scale=[];
        scale2=[];
        % figure
        
        Xstart=1;
        for c=1:length(cond)
            temp=[];
            
            
            subplot(3,6,1:5)
            hold on
            temp=nan(strides(c),size(X1,2));
            %         if c==3
            %          temp(1:sz(c),:)=movmean(X1(strides(c)-40:strides(c),:),5,'omitnan'); %moving average of 5
            %         else
            
            temp(1:sz(c),:)=movmean(X1(ini(c):strides(c),:),5,'omitnan'); %moving average of 5
            %         end
            %         y=[];
            y=nanmean(temp,2)'; %across subjects mean
            %         y(isnan(y))=[];
            
            if max(abs(y))>1  %To scale the hidden state to be max 1
                
                scale=1;%(1/max(abs(y)));
                
            end
            
            if isempty(scale)
                scale=1;
            end
            condLength=length(y);
            if c==4
                if group==1 || group==2 || group==3 
                    x=1010:1010+condLength-1;
                else
                    x=680:680+condLength-1;
                end
            else
                x=Xstart:Xstart+condLength-1;
            end
            E=std(temp,0,2,'omitnan')./sqrt(size(temp,2)); %Standar error
            %         E(isnan(E))=[];
            [Pa, Li{loop}]= nanJackKnife(x,y*scale,E',color,color+0.5.*abs(color-1),Opacity);
            yline(0)
            
            if c==4
                subplot(3,6,6)
                bar(loop,y(1)*scale,'FaceColor',color,'BarWidth',0.9)
                hold on
                errorbar(loop,y(1)*scale,E(1),'LineStyle','none','color','k','LineWidth',2)
                for ss=1:size(X1,2)
                    plot(loop-0.03,temp(1,ss)*scale,'.','MarkerSize', 15, 'Color',[150 150 150]./255)
                end
                
            end
            
            %         bar(plotHandles(p),xval(:,i),squeeze(plotData(i,p,:)),'FaceColor',colors(i,:),'BarWidth',0.2)
            %         errorbar(plotHandles(p),xval(:,i),squeeze(plotData(i,p,:)),squeeze(varData(i,p,:)),'LineStyle','none','color','k','LineWidth',2)
            
            
            subplot(3,6,7:11)
            hold on
            temp=nan(strides(c),size(X2,2)); %moving average of 5
            temp(1:sz(c),:)=movmean(X2(ini(c):strides(c),:),5,'omitnan'); %across subjects mean
            y=nanmean(temp,2)';
            %        y(isnan(y))=[];
            
            if max(abs(y))>1  %To scale the hidden state to be max 1
                
                scale2=1;%(1/max(abs(y)))
                
            end
            
            if isempty(scale2)
                scale2=1;
            end
            
            condLength=length(y);
            
            %         x=Xstart:Xstart+condLength-1;
            E=std(temp,0,2,'omitnan')./sqrt(size(temp,2)); %Standar error
            %         E(isnan(E))=[];
            [Pa, Li{loop}]= nanJackKnife(x,y*scale2,E',color,color+0.5.*abs(color-1),Opacity);
            yline(0)
            
            if c==4
                subplot(3,6,12)
                bar(loop,y(1)*scale2,'FaceColor',color,'BarWidth',.9)
                hold on
                errorbar(loop,y(1)*scale2,E(1),'LineStyle','none','color','k','LineWidth',2)
                
                for ss=1:size(X1,2)
                    plot(loop-0.03,temp(1,ss)*scale2,'.','MarkerSize', 15, 'Color',[150 150 150]./255)
                end
            end
            
            subplot(3,6,13:17)
            hold on
            temp=nan(strides(c),size(X2,2)); %moving average of 5
            temp(1:sz(c),:)=movmean(X3(ini(c):strides(c),:),5,'omitnan'); %across subjects mean
            y=nanmean(temp,2)';
            %        y(isnan(y))=[];
            
            if max(abs(y))>1  %To scale the hidden state to be max 1
                
                scale2=1;%(1/max(abs(y)));
                
            end
            
            if isempty(scale2)
                scale2=1;
            end
            
            condLength=length(y);
            
            %         x=Xstart:Xstart+condLength-1;
            E=std(temp,0,2,'omitnan')./sqrt(size(temp,2)); %Standar error
            %         E(isnan(E))=[];
            [Pa, Li{loop}]= nanJackKnife(x,y*scale2,E',color,color+0.5.*abs(color-1),Opacity);
            yline(0)
            
            if c==4
                subplot(3,6,18)
                bar(loop,y(1)*scale2,'FaceColor',color,'BarWidth',.9)
                hold on
                errorbar(loop,y(1)*scale2,E(1),'LineStyle','none','color','k','LineWidth',2)
                
                for ss=1:size(X1,2)
                    plot(loop-0.03,temp(1,ss)*scale2,'.','MarkerSize', 15, 'Color',[150 150 150]./255)
                end
            end
            
            
            Xstart=Xstart+condLength;
            
            hold on
            
            
        end
        %%
        
        subplot(3,6,1:5)
        if group==1 || group==2 || group==3
            if length(cond)==4
                pp=patch([50 1000 1000 50],[-1.2 -1.2 1.2 1.2],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
                uistack(pp,'bottom')
            else
                
                pp=patch([50 950 950 50],[-1.2 -1.2 1.2 1.2],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
                uistack(pp,'bottom')
            end
            ylabel('X_{reactive}')
            axis tight
        else
            if length(cond)==4
                pp=patch([50 650 650 50],[-1.2 -1.2 1.2 1.2],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
                uistack(pp,'bottom')
            else
                pp=patch([50 600 600 50],[-1.2 -1.2 1.2 1.2],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
                uistack(pp,'bottom')
            end
            
            ylabel('X_{reactive}')
            axis tight
            
        end
        
        
        subplot(3,6,7:11)
        if group==1 || group==2 || group==3 
            
            if length(cond)==4
                pp=patch([50 1000 1000 50],[-1.2 -1.2 1.2 1.2],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
                uistack(pp,'bottom')
            else
                
                pp=patch([50 950 950 50],[-1.2 -1.2 1.2 1.2],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
                uistack(pp,'bottom')
            end
            
            ylabel('X_{reactive}')
            axis tight
        else
            if length(cond)==4
                pp=patch([50 650 650 50],[-1.2 -1.2 1.2 1.2],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
                uistack(pp,'bottom')
            else
                pp=patch([50 600 600 50],[-1.2 -1.2 1.2 1.2],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
                uistack(pp,'bottom')
            end
            
            ylabel('X_{reactive}')
            axis tight
            
        end
        
        
        subplot(3,6,13:17)
        if group==1 || group==2 || group==3
            
            if length(cond)==4
                pp=patch([50 1000 1000 50],[-1.2 -1.2 1.2 1.2],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
                uistack(pp,'bottom')
            else
                
                pp=patch([50 950 950 50],[-1.2 -1.2 1.2 1.2],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
                uistack(pp,'bottom')
            end
            
            ylabel('X_{reactive}')
            axis tight
        else
            if length(cond)==4
                pp=patch([50 650 650 50],[-1.2 -1.2 1.2 1.2],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
                uistack(pp,'bottom')
            else
                pp=patch([50 600 600 50],[-1.2 -1.2 1.2 1.2],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
                uistack(pp,'bottom')
            end
            
            ylabel('X_{reactive}')
            axis tight
            
        end
        
        subplot(3,6,6)
        set(gca, 'XTick', 1:length(groups), 'XTickLabels', labels(groups))
        subplot(3,6,12)
        set(gca, 'XTick', 1:length(groups), 'XTickLabels', labels(groups))
        
        set(gcf,'color','w')
        legend([Li{:}]',[labels(groups)])
    end
    
end