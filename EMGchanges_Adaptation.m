
function EMGchanges_Adaptation(normalizedGroupData,newLabelPrefix,ep,fi,ph,baseEp,fdr,summFlag,groupID,savefig)

% Code to plot the changes in EMG during early and late adaptation 


for ee=1:length(ep)

[fi,pc,labels,dataEc,dataRefc]=plotCheckerboards(normalizedGroupData,newLabelPrefix,ep(ee,:),fi,ph(1,ee),baseEp,2,summFlag);

%%
%nonparametric stats
dataEcmed=transpose(squeeze(nanmedian(dataEc,4)));
%these values are generated to assess if effect sizes are larger than threshold value
%-this matrix is in the same format as the checkerboards (rows are muscles and columns are gait phases
%-note that the first row corresponds to the bottom row of the
% checkerboards given that x-values in plot go from 1:30

[pvalc,hc,alphaAdj_c]=checkerstatsV2(dataEc,[],1,0,fdr,'benhoch',0);%mindif has to be zero, since signrank cannot reliably do a two-tail test agains another value


%-matrices hc and pvalc are in the same format as the checkerboards (see dataEcmed)
disp(['p-threshold controls = ',num2str(alphaAdj_c)])
get(pc);
hold on
for i=1:size(hc,1)
    for k=1:size(hc,2)
        if hc(i,k)==1  && abs(dataEcmed(i,k))>0.2 %since statistical testing was done againts zero, amplitude testing happens here 10% change from baseline
            plot3((k-0.5)/12,i-0.5,1,'.','MarkerSize',15,'Color','k');
           
        end
    end
end
title([cell2mat(ep.Properties.ObsNames(ee)),' p=',num2str(round(1000*alphaAdj_c)/1000)]);
end


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
    0.6350      0.0780    0.1840];
ex1=cc(2,:);
ex2=cc(5,:);
mid=ones(1,3);
N=100;
gamma=1.5; %gamma > 1 expands the white (mid) part of the map, 'hiding' low values. Gamma<1 does the opposite
gamma=1;
map=[flipud(mid+ (ex1-mid).*([1:N]'/N).^gamma); mid; (mid+ (ex2-mid).*([1:N]'/N).^gamma)];


%%
fs=14;
colormap(flipud(map))
% colormap default
set(gcf,'color','w');                                     
set(ph(:,1),'CLim',[-1 1]*1,'FontSize',fs);
set(ph(:,2:end),'YTickLabels',{},'CLim',[-1 1]*1,'FontSize',fs);
colorbar 
set(gcf,'renderer','painters')
if savefig==1
    saveas(gcf,[groupID, '_EMGstructure.eps'])
end
end
%%

