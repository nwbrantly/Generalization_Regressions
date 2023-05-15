%creating subspace 

load BATR_5_AsymC13_ShortPertubations_RemovedBadMuscle_0.mat
fname= 'dynamicsData_BATR_subj_5_RemoveBadMuscles0_splits_0.h5';



Casym=[ C(:,5)/norm(C(:,5)) C(:,6)/norm(C(:,6))];% C(:,10)/norm(C(:,10))];
r =rand(10000,size(Casym,2));
Subspace= Casym *r';

EMGdata=h5read(fname,'/EMGdata');
[Y,Yasym,~,U,~,Ysum]=groupDataToMatrixForm(1:size(EMGdata,3),0,fname);

Uf=[U;ones(size(U))];
removebaseline=1
if removebaseline==1
    bias=nanmean(Yasym(5:30,:));
    C=Casym-bias';
    Ymodel=Yasym'-bias'; % Transpose the EMG data to match equations
end
[pp,cc,aa]=pca(Subspace','Centered','off');
% [coeff,score,latent,tsquared,explained,mu]=pca(Yasym,'Centered','off');
%%Input has to be row observation  and columns variables 
%%pp - data projected in the principal component 
%%cc - Corresponding matrix of eigenvectors
%%aa - vector of eigent values 

figure 
subplot(2,1,1)
plot(aa,'LineWidth',2)
ylabel('Eigenvalues')
xline(2,'r')
xlabel('Number of Components')
% ylim([-.5 7])
 axis tight

% yline(0.85,'r')
% xline(0.85,'r')

subplot(2,1,2)
csum=cumsum(aa);
variance_explained = csum / sum(aa);
plot(variance_explained,'LineWidth',2)
xline(2,'r')
ylabel('Variance explained')
xlabel('Number of Components')
set(gcf,'color','w')

Cres=[pp(:,1:2)];

% Cres=[pp(:,1)+pp(:,3) pp(:,2) pp(:,4) ];
Cinv=pinv(Cres);
model.C=Cres;
Xresidual= Yasym*Cinv';



binwith=5;
figure
scatter(1:length(movmean(Xresidual(:,1),binwith)),movmean(Xresidual(:,1),binwith),10,'k','filled')
hold on 
scatter(1:length(movmean(Xresidual(:,2),binwith)),movmean(Xresidual(:,2),binwith),10,'b','filled')

model.C=Cres;
legacy_vizSingleModel_FreeModel_ShortAdaptation(model,Yasym',Uf,0)
%%

figure
plot3(Casym(:,1),Casym(:,2),Subspace)
%%

figure
hold on
for s=1:200;%size(r,1)
    
    scatter3(Casym(:,1),Casym(:,2),Subspace(:,s),'b','filled')
    
end


%%

for s=1:680
temp(s)=norm(Yasym(s,:));
end