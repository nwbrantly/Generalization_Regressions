%  load('ATS_fixD_ShortPostAdapt120422T194824.mat')
 close all
 Y=datSet.out;
%  Y=Y(:,1:1150);
% [pp,cc,aa]=pca(Y');
% [pp,cc,aa]=pca(Y');
[pp,cc,aa]=pca(Y','Centered','off');
%%pp - data projected in the principal component 
%%cc - Corresponding matrix of eigenvectors
%%aa - vector of eigent values 

%Recontruction of the data 
Xcentered = pp*cc'; 
Ymean=nanmean(Y,2);
Ynew=Xcentered;%+Ymean;

figure 
subplot(2,1,1)
plot(aa,'LineWidth',2)
ylabel('Eigenvalues')
xline(3,'r')
xlabel('Number of Components')
% ylim([-.5 7])
 axis tight
% yline(0.85,'r')
% xline(0.85,'r')

subplot(2,1,2)
csum=cumsum(aa);
variance_explained = csum / sum(aa);
plot(variance_explained,'LineWidth',2)
xline(3,'r')
ylabel('Variance explained')
xlabel('Number of Components')
set(gcf,'color','w')
%% Using first 2 components to estimate the states 

C=[pp(:,1:3)];
Cinv=pinv(C)';
X= Y'*Cinv;

figure
subplot(3,1,1)
plot(movmean(X(:,1),5),'LineWidth',2)
yline(0)
title('C_1 = PCA1')
% legend('ATR','ATS')

subplot(3,1,2)

plot(movmean(X(:,2),5),'LineWidth',2)
yline(0)
title('C_2 = PCA2')

subplot(3,1,3)
plot(movmean(X(:,3),5),'LineWidth',2)
yline(0)
title('C_3 = PCA3')
% legend('ATR','ATS')

set(gcf,'color','w')
%%
model.C=C;
legacy_vizSingleModelMLMC_FreeModel(model,datSet.out,datSet.in)

model.C=C(:,1:2);
legacy_vizSingleModelMLMC_FreeModel(model,datSet.out,datSet.in)