%colors
    poster_colors;
    colorOrder=[[0.4940 0.1840 0.5560];[0.9290 0.6940 0.1250]; [0.4660 0.6740 0.1880]];
    
    
   
%     figure

for t=1:2
    figure
    for g=1:2
        
        if g==1
            load('/Users/dulcemariscal/Documents/GitHub/Generalization_Regressions/RegressionAnalysis/RegModelResults_ATR01-03_V2/GroupResults/ATRdefaultsplit_1flip_1_group_models_ver00.mat')
        elseif g==2
            load('/Users/dulcemariscal/Documents/GitHub/Generalization_Regressions/RegressionAnalysis/RegModelResults_ATRS01_V2/ATS01defaultsplit_1models_ver00.mat')
            
        end
         transition={fitTrans1NoConst, fitTrans2NoConst};
        for b=1:3
            if g==1
                x=[1,2,3];
            else
                x=[5,6,7];
            end
            hold on
            if g==1
             bar(x(b),abs(transition{t}.Coefficients.Estimate(b)),'FaceColor',colorOrder(b,:))
            else
                bar(x(b),abs(transition{t}.Coefficients.Estimate(b)),'FaceColor','w','EdgeColor',colorOrder(b,:))
            end
            errorbar(x(b),abs(transition{t}.Coefficients.Estimate(b)),fitTrans1NoConst.Coefficients.SE(b),'k')
            
        end
        title(['Transition', num2str(t)])
        ylabel('|\beta|')
        xticks([1 2 3 5 6 7])
        xticklabels({'\beta_{adapt}','\beta_{non-adapt}','\beta_{env}','\beta_{adapt}','\beta_{non-adapt}','\beta_{env}'})
        set(gcf,'color','w');
        legend('Learning','Generalization')
        axis([0 8 -0.1 1.5])
    end
end

%%

groupID={'NTS','NTR'};
colorOrder=[[0.4940 0.1840 0.5560];[0.9290 0.6940 0.1250]; [0.4660 0.6740 0.1880]];

for i=1:2
    load(['GroupRegression00_',groupID{i},'.mat'])
    
    
    figure
    
    for tr=1:3
        hold on
        eval(['data= trans',num2str(tr)]);
        
        betas_order={1:3, 5:7, 9:11};
        c=0;
        for beta=betas_order{tr}
            c=c+1;
            if (contains(groupID{i}, 'TR'))
                b{c}=bar(beta,abs(data.Coefficients.Estimate(c)), 'FaceColor' ,colorOrder(c,:));
             
            elseif (contains(groupID{i}, 'TS'))
                
                b{c}=bar(beta,abs(data.Coefficients.Estimate(c)), 'FaceColor' ,'w','EdgeColor',colorOrder(c,:));
%                 errorbar(beta, abs(data.Coefficients.Estimate(c)),abs(trans1.Coefficients.SE(c)),'k')
            end
               errorbar(beta, abs(data.Coefficients.Estimate(c)),abs(trans1.Coefficients.SE(c)),'k')
            
        end
        
    end
    
    title(groupID{i})
    legend([b{:}],'\beta_{adapt}','\beta_{noadapt}','\beta_{env}')
    ylabel('|\beta|')
    xticks([2 6 10])
    xticklabels({'Trans 1: Long Split','Trans 2: Env Changes','Trans 3: Short Split '})
    set(gcf,'color','w');
    
    
    
end
