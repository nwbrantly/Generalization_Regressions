                                                                                                              %Transition 1 
%% R^2 ordinary 
 load('GroupRegressionParams.mat')  
 figure
 hold on
 plot(1,Rsquared1_Ord(1),'or','MarkerFaceColor','r')
 
 load('IndivRegressionParams.mat') 
  plot(1.05,Rsquared1_Ord(:,1),'ok')
 
 load('GroupRegressionParamsYA.mat')
 plot(1.5,Rsquared1_Ord(1),'om','MarkerFaceColor','m')
 load('IndivRegressionParamsYA.mat')
 plot(1.55,Rsquared1_Ord(:,1),'om')
 
ylabel('R^2, Ordinary')
title('Transition 1')
xticks([1 1.5])
xticklabels({'OA','ATR'})
axis([0.9 1.8 0 0.8])
set(gcf,'color','w');
%% R^2 adjusted
 load('GroupRegressionParams.mat')  
 figure
 hold on
 plot(1,Rsquared1_Adj(1),'or','MarkerFaceColor','r')
 
 load('IndivRegressionParams.mat') 
  plot(1.05,Rsquared1_Adj(:,1),'ok')
 
 load('GroupRegressionParamsYA.mat')
  plot(1.5,Rsquared1_Adj(1),'om','MarkerFaceColor','m')
 load('IndivRegressionParamsYA.mat')
 plot(1.55,Rsquared1_Adj(:,1),'om')
 
ylabel('R^2, Adjusted')
title('Transition 1')
xticks([1 1.5])
xticklabels({'OA','ATR'})
axis([0.9 1.8 0 0.8])
set(gcf,'color','w');
%% Betas 

 load('GroupRegressionParams.mat')  
 figure
 hold on
 plot(1,Coefficients_trans1(1),'or','MarkerFaceColor','r')
 plot(2,Coefficients_trans1(2),'or','MarkerFaceColor','r')
 plot(3,Coefficients_trans1(3),'or','MarkerFaceColor','r')
 
 load('IndivRegressionParams.mat') 
 plot(1.05,Coefficients_trans1(1,:),'ok')
 plot(2.05,Coefficients_trans1(2,:),'ok')
 plot(3.05,Coefficients_trans1(3,:),'ok')
 
 load('GroupRegressionParamsYA.mat')
 plot(1.1,Coefficients_trans1(1),'om','MarkerFaceColor','m')
 plot(2.1,Coefficients_trans1(2),'om','MarkerFaceColor','m')
 plot(3.1,Coefficients_trans1(3),'om','MarkerFaceColor','m')
  
  
 load('IndivRegressionParamsYA.mat')
 plot(1.15,Coefficients_trans1(1,:),'om')
 plot(2.15,Coefficients_trans1(2,:),'om')
 plot(3.15,Coefficients_trans1(3,:),'om')
 
ylabel('\beta')
title('Transition 1')
xticks([1 2 3])
xticklabels({'\beta_{adapt}','\beta_{non-adapt}','\beta_{env-switch}'})
axis([0.5 3.5 -1 1.5])
set(gcf,'color','w');
%%
%Transition 2
%% R^2 ordinary 
 load('GroupRegressionParams.mat')  
 figure
 hold on
 plot(1,Rsquared2_Ord(1),'or','MarkerFaceColor','r')
 
 load('IndivRegressionParams.mat') 
  plot(1.05,Rsquared2_Ord(:,1),'ok')
 
 load('GroupRegressionParamsYA.mat')
 plot(1.5,Rsquared2_Ord(1),'om','MarkerFaceColor','m')
 load('IndivRegressionParamsYA.mat')
 plot(1.55,Rsquared2_Ord(:,1),'om')
 
ylabel('R^2, Ordinary')
title('Transition 2')
xticks([1 1.5])
xticklabels({'OA','ATR'})
axis([0.9 1.8 0 0.8])
set(gcf,'color','w');
%% R^2 adjusted
 load('GroupRegressionParams.mat')  
 figure
 hold on
 plot(1,Rsquared2_Adj(1),'or','MarkerFaceColor','r')
 
 load('IndivRegressionParams.mat') 
  plot(1.05,Rsquared2_Adj(:,1),'ok')
 
 load('GroupRegressionParamsYA.mat')
  plot(1.5,Rsquared2_Adj(1),'om','MarkerFaceColor','m')
 load('IndivRegressionParamsYA.mat')
 plot(1.55,Rsquared2_Adj(:,1),'om')
 
ylabel('R^2, Adjusted')
title('Transition 2')
xticks([1 1.5])
xticklabels({'OA','ATR'})
axis([0.9 1.8 0 0.8])
set(gcf,'color','w');
%% Betas 

 load('GroupRegressionParams.mat')  
 figure
 hold on
 plot(1,Coefficients_trans2(1),'or','MarkerFaceColor','r')
 plot(2,Coefficients_trans2(2),'or','MarkerFaceColor','r')
 plot(3,Coefficients_trans2(3),'or','MarkerFaceColor','r')
 
 load('IndivRegressionParams.mat') 
 plot(1.05,Coefficients_trans2(1,:),'ok')
 plot(2.05,Coefficients_trans2(2,:),'ok')
 plot(3.05,Coefficients_trans2(3,:),'ok')
 
 load('GroupRegressionParamsYA.mat')
 plot(1.1,Coefficients_trans2(1),'om','MarkerFaceColor','m')
 plot(2.1,Coefficients_trans2(2),'om','MarkerFaceColor','m')
 plot(3.1,Coefficients_trans2(3),'om','MarkerFaceColor','m')
  
  
 load('IndivRegressionParamsYA.mat')
 plot(1.15,Coefficients_trans2(1,:),'om')
 plot(2.15,Coefficients_trans2(2,:),'om')
 plot(3.15,Coefficients_trans2(3,:),'om')
 
ylabel('\beta')
title('Transition 2')
xticks([1 2 3])
xticklabels({'\beta_{adapt}','\beta_{non-adapt}','\beta_{env-switch}'})
axis([0.5 3.5 -1 1.5])
set(gcf,'color','w');