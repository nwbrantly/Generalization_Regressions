function [fh,model] = vizSingleModel_LinearRegression_IndvLeg(singleModel,Y,lower,labels,isf)

%Funtion to plot the time courses of the hidden states, data, fit and the
%residual of the model
%This function was initially created to plot 


%also see legacy_vizSingleModelMLMC

if isa(singleModel,'linsys')
    singleModel=struct(singleModel);
end

M=size(zeros(2,2),1);
fh=figure('Units','Normalized','OuterPosition',[0 0 .5 1],'Color',ones(1,3));
if nargin>1
    Nx=12;
else
    Nx=4;
end
singleModel.B=zeros(2,2);
Nu=size(singleModel.B,2);
% Nu=[2];
Ny=M+1+Nu;
model{1}=singleModel;
Nu=size(model{1}.B,2);
Nc=size(Y,1);

if nargin<3
    lower=false;
end
%% First: Define which model you want to run 
if lower==1 % Using reactive and adaptation without dynamics 
    
    C=model{1}.C;
    if size(C,2)==1
         X=1*(ones(size(C,2),size(Y,2)))' ;
    else
        
        X=0.5*(ones(size(C,2),size(Y,2)))' ;
    end

    Y2= C * X' ; % Data reconstructed with the perdetermine dynamics
    model{1}.Data=Y;
    model{1}.States=X;
    model{1}.Out=Y2;
    model{1}.Res=Y-Y2;
    

elseif lower==2 %Analysis using PCA in all the data 
    [pp,cc,aa]=pca(Y');
    NNMF= [(cc(:,1:2)*pp(:,1:2)') + nanmean(Y')]';
    X=cc(:,1:2);
    C=pp(:,1:2);

    model{1}.Data=Y;
    model{1}.States=cc(:,1:2);
    model{1}.Out=NNMF;
    model{1}.Res=Y - NNMF;
    
elseif lower==3 %removing the mean on the data. similar analysis as PCA 
    m=nanmean(Y,2); %m stands fro mean
    model{1}.C=model{1}.C-m;
    Ydmean=Y-m;
    
    C=model{1}.C; %getting C from the struture
    Cinv=pinv(C)'; %pseudoinverse of the C (C is not a squared matrix)
    X = Ydmean'*Cinv;
    Y2=  model{1}.C * X' ; %yhat = C 
    
    model{1}.States = X; %x= y/C
    model{1}.Out=   Y2 + m ; %yhat = C
    model{1}.Res=Y-model{1}.Out;
    model{1}.Data=Y;
    
else
    C=model{1}.C; %getting C from the struture
    Cinv=pinv(C)'; %pseudoinverse of the C (C is not a squared matrix)
    X = Y'*Cinv;  %x= y/C getting the dynamics of the hidden states (we are using least-sqaures seee: Penrose, Roger (1956))
    Y2= C * X' ;
    model{1}.Data=Y;
    model{1}.States=X;
    model{1}.Out=Y2;
    model{1}.Res=Y-Y2;
end
    
%% Define colormap:
% ex1=[1,0,0];
% ex2=[0,0,1];
% cc=[0    0.4470    0.7410
%     0.8500    0.3250    0.0980
%     0.9290    0.6940    0.1250
%     0.4940    0.1840    0.5560
%     0.4660    0.6740    0.1880
%     0.3010    0.7450    0.9330
%     0.6350    0.0780    0.1840];
% ex1=cc(2,:);
% ex2=cc(5,:);
ex2=[0.2314    0.2980    0.7529];
ex1=[0.7255    0.0863    0.1608];
mid=ones(1,3);
N=100;
gamma=1.5; %gamma > 1 expands the white (mid) part of the map, 'hiding' low values. Gamma<1 does the opposite
gamma=1;
map=[flipud(mid+ (ex1-mid).*([1:N]'/N).^gamma); mid; (mid+ (ex2-mid).*([1:N]'/N).^gamma)];

%% Plot dynamics and heatmaps
% ytl={'HIP','GLU','TFL','RF','VL','VM','SEMT','SEMB','BF','MG','LG','SOL','PER','TA'};
%
if ~isempty(labels)
  muscleOrder=  labels;
else
muscleOrder={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT', 'VM', 'VL', 'RF', 'HIP','TFL', 'GLU'};
end

ytl=([strcat('f',muscleOrder) strcat('s',muscleOrder)]);  %List of muscle 

ytl(end:-1:1) = ytl(:);
if isf==1
ytl=ytl(length(ytl(:))\2+1:end); %fast
end

yt=1:length(ytl);
fs=7;
% STATES
C=[model{1}.C];
X=[X'];


rotMed='none'; % I am not rotating the data 
% [CDrot,XUrot]=rotateFac(CD,XU,rotMed);
Crot=C;
XUrot=X;

if strcmp(rotMed,'none')
%     factorName=[strcat('C_',num2str([1:size(model{1}.C,2)]'));strcat('D_',num2str([1:size(model{1}.D,2)]'))];
    factorName=[strcat('C_',num2str([1:size(model{1}.C,2)]'))];
%     latentName=[strcat('State ',' ', num2str([1:size(model{1}.C,2)]'));strcat('Input ',' ', num2str([1:size(model{1}.D,2)]'))];
    latentName=[strcat('State ',' ', num2str([1:size(model{1}.C,2)]'))];
else
    factorName=strcat('Factor ',num2str([1:size(C,2)]'));
    latentName=strcat('Latent ',num2str([1:size(C,2)]'));
end

aC=prctile(abs(Y(:)),98);

projY=X;
for i=1:size(C,2)
   
    hold on
    if i<=size(C,2)
        ph(i)=subplot(Nx,size(C,2),i); %TOP row: states temporal evolution and data projection
        scatter(1:size(Y,2),projY(i,:),5,.7*ones(1,3),'filled')
        title(latentName(i,:))
        p(i)=plot(XUrot(i,:),'LineWidth',2,'Color','k');
        ax=gca;
        ax.Position=ax.Position+[0 .045 0 0];
        axis tight
        grid on
        set(gca,'ColorOrderIndex',1)
    end
   
    subplot(Nx,Ny,Ny+i+[0,Ny])% Second row: checkerboards
    try
        imagesc((reshape(Crot(:,i),12,Nc/12)'))
    catch
        imagesc((Crot(:,i)))
    end
    set(gca,'XTick',[],'YTick',yt,'YTickLabel',ytl,'FontSize',fs)
    ax=gca;
    ax.YAxis.Label.FontSize=12;
    colormap(flipud(map))
%     colormap(map)
%     caxis([-aC aC])
    caxis([-1 1])
    axis tight
    title(factorName(i,:))
    ax=gca;
    ax.Position=ax.Position+[0 .03 0 0];
    hold on
    aa=axis;
    plot([0.1 1.9]+.5, [1.8 1.8],'k','LineWidth',2,'Clipping','off')
    text(.8,1.6,'DS','FontSize',6)
    plot([2.1 5.9]+.5, [1.55 1.55],'k','LineWidth',2,'Clipping','off')
    text(2.8,1.6,'SINGLE','FontSize',6)
    plot([6.1 7.9]+.5, [1.55 1.55],'k','LineWidth',2,'Clipping','off')
    text(6.9,1.6,'DS','FontSize',6)
    plot([8.1 11.9]+.5, [1.55 1.55],'k','LineWidth',2,'Clipping','off')
    text(9,1.6,'SWING','FontSize',6)
    axis(aa)
    

end
subplot(Nx,Ny,Ny+i+1+[0,Ny])% This one is display just for the colorbar 
caxis([-1 1])
colorbar
colorbar('Ticks',[0,1],'TickLabels',{'0','100%'})
ax=gca;
ax.Position=ax.Position+[0 .03 0 0];
hold on
aa=axis;
axis(aa)

% linkaxes(ph,'y')
%Covariances
%subplot(Nx,Ny,Ny)
%imagesc(model{1}.Q)
%set(gca,'XTick',[],'YTick',[],'YTickLabel',ytl,'FontSize',8)
%colormap(flipud(map))
%aC=.5*max(abs(model{1}.Q(:)));
%caxis([-aC aC])
%axis tight
%subplot(Nx,Ny,2*Ny+[0,Ny])
%imagesc(model{1}.R)
%set(gca,'XTick',[],'YTick',[],'YTickLabel',ytl,'FontSize',8)
%colormap(flipud(map))
%aC=.5*max(abs(model{1}.R(:)));
%caxis([-aC aC])
%axis tight

if nargin<2
    %Third row: one-ahead step-response


else %IF DATA PRESENT:
N=size(Y,2);
viewPoints=[3,35,44,450,483,675,682]; %PATR - PATS 
% viewPoints=[3,437,443,635]; %PATR - PATS 
binw=4; %Plus minus 2
viewPoints(viewPoints>N-binw/2)=[];
Ny=length(viewPoints);
M=length(model);
% dd=Y(:,2:39);
% dd=dd-mean(dd,2); %Baseline residuals under flat model
meanVar=1;%mean(sum(dd.^2,1),2);
for k=1:3
    for i=1:Ny
        switch k
        case 1 % Third row, actual data
            dd=Y(:,viewPoints(i)+[-(binw/2):(binw/2)]);
            nn='data';
        case 2 %Fourth row: one-ahead data predictions
            dd=model{1}.Out(:,viewPoints(i)+[-(binw/2):(binw/2)]);
            nn={'Data Fit'};
        case 3 % Fifth row:  data residuals (checkerboards)
            dd=model{1}.Res(:,viewPoints(i)+[-(binw/2):(binw/2)]);
            nn='residual';
        end

        subplot(Nx,Ny,i+(1+2*k)*Ny+[0,Ny])
        try
            imagesc(reshape(median(dd,2),12,size(Y,1)/12)')

        catch
            imagesc(median(dd,2))
        end
        
        set(gca,'XTick',[],'YTick',yt,'YTickLabel',ytl,'FontSize',fs)

        ax=gca;

        colormap(flipud(map))
%         caxis([-aC aC])
        caxis([-1 1])
        axis tight
        if k==1
            title(['Output at t=' num2str(viewPoints(i))])
%             txt={'Base Late','Early Adapt','Late Adap','Early Post','Mid Post','Late Post'};
%             title(txt{i})
            ax=gca;
            ax.Title.FontSize=10;
        end
        if k==2
            ax=gca;
            ax.Title.FontSize=10;
        end
        
        
        if k==3
%             title(['RMSE=' num2str(sqrt(mean(sum(dd.^2,1),2)))])
            title(['RMSE=' num2str(sqrt(mean(mean((dd.^2),1),2))/sqrt(meanVar))]);
%              title(['normalized RMSE=' num2str(sqrt(mean(sum(dd.^2,1),2))/sqrt(meanVar))])
%              title(['RRMSE=' num2str(sqrt(mean((dd.^2),1)))])
%             title(['RMSE=' num2str(sqrt(mean(mean(dd,2).^2)))])
        end
        if i==1
            ylabel(nn)
            ax=gca;
            ax.YAxis.Label.FontWeight='normal';
            ax.YAxis.Label.FontSize=12;
        end
    end
end


% PCA analysis for upper bound 
% Y2=Y;
% adapt=Y(:,41:480);
% post=Y(:,481:491);
% [pp,cc,aa]=pca(adapt','Centered','off');
% [pp_2,cc_,aa_]=pca(post','Centered','off');
% [coeff,score,latent,tsquared,explained,mu]=pca(Yasym,'Centered','off');
%%Input has to be row observation  and columns variables 
%%pp - Corresponding matrix of eigenvectors
%%cc - data projected in the principal component 
%%aa - vector of eigent values
% 
% W=[pp(:,1:2) pp_2(:,1)];
% Winv=pinv(W);
% dynamics= Y'*Winv';
% PC= W * dynamics' ; 
% Res= Y - PC;
% 
% [pp,cc,aa]=pca(Y');%,'Centered','off');
% PC= [(cc(:,1:2)*pp(:,1:2)') + nanmean(Y')]'; 
% [pp,cc,aa]=pca(Y','Centered','off');
% Winv=pinv(W);
% dynamics= Y'*Winv';
% PC= W * dynamics' ;

% NNMF
% YA=[model{1}.Data(:,1:70) model{1}.Data(:,450:481+60)];
% YA=[model{1}.Data(:,1:70) model{1}.Data(:,450:480)];
% % swift=abs(min(model{1}.Data',[],'all'));
% % data=YA+swift;
% data=YA;
% % data2=model{1}.Data+swift;
% data2=model{1}.Data;
% [W,H] = nnmf(data,4);
% Casym = W';
% Cinv=pinv(Casym); %Gettign the pseudoinverse of the EMGreactive and EMGcontext
% % Wasym_NNMF = Cinv*Y'; %x= y/C
% Wasym_NNMF_noshift = Cinv'*Y; %x= y/C
% NNMF=  (Casym'* Wasym_NNMF_noshift); %yhat = C 
% ResNNMF= Y- NNMF;

% Res=PC- movmean(Y,k);

% 2 states 
% C=C(:,1:2); %getting C from the struture
% Cinv=pinv(C)'; %pseudoinverse of the C (C is not a squared matrix)
% X = Y'*Cinv; %x= y/C getting the dynamics of the hidden states (we are using least-sqaures seee: Penrose, Roger (1956))
% Y_2= C* X' ; % Data reconstructed with the perdetermine dynamics
% Res_2=Y-Y_2;
% 


%% RMSE and R2
%% RMES plot
Ny=1;
subplot(Nx,Ny,1+9*Ny)
hold on
binw=5;

dd=model{1}.Res;
aux1=sqrt(mean(dd.^2))/sqrt(meanVar);
aux1=conv(aux1,ones(1,binw)/binw,'valid'); %Smoothing
plot(aux1,'LineWidth',2,'DisplayName',[num2str(size(C,2)), 'states'],'Color','r');
 
% aux2=sqrt(sum(Res.^2))/sqrt(meanVar);
% aux2=conv(aux2,ones(1,binw)/binw,'valid'); %Smoothing
% plot(aux2,'LineWidth',2,'DisplayName','PCA ','Color','k');

% aux2=sqrt(mean(ResNNMF.^2))/sqrt(meanVar);
% aux2=conv(aux2,ones(1,binw)/binw,'valid'); %Smoothing
% plot(aux2,'LineWidth',2,'DisplayName','NNMF','Color','k');%,"#A2142F");
% 
% aux3=sqrt(sum(Res_2.^2))/sqrt(meanVar);
% aux3=conv(aux3,ones(1,binw)/binw,'valid'); %Smoothing
% plot(aux3,'-.','LineWidth',2,'DisplayName','2 states','Color','b');
% 
% % Instantaneous SD
Res_3= conv2(Y,[0,1,-1],'valid'); %Y(k)-.5*(y(k+1)+y(k-1));
% Res_3= conv2(Yasym,[0,1,-1],'valid'); %Y(k)-.5*(y(k+1)+y(k-1));
aux3=sqrt(mean(Res_3.^2))/sqrt(2);
aux3=conv(aux3,ones(1,binw)/binw,'valid'); %Smoothing
plot([nan nan aux3],'-','LineWidth',2,'DisplayName','Instantaneous SD','Color',"#0072BD");
ylabel({'residual';' RMSE'})

% aux3=sqrt(sum(Res_4.^2))/sqrt(k);
% aux3=conv(aux3,ones(1,binw)/binw,'valid'); %Smoothing
% plot(aux3,'-.','LineWidth',2,'DisplayName',['2 STD'],'Color',"#D95319");
% ylabel({'residual';' RMSE'})
%title('MLE one-ahead output error (RMSE, mov. avg.)')

%Add previous stride model:
% ind=find(diff(U(1,:))~=0);
% % Y(:,ind)=nan;
% aux1=(Y(:,2:end)-Y(:,1:end-1))/sqrt(2);
% aux1=sqrt(mean(aux1.^21sqrt(meanVar);
% % aux1=sqrt(mean(aux1.^2));%/sqrt(meanVar);
% aux1=conv(aux1,ones(1,binw)/binw,'valid'); %Smoothing
% % plot(aux1,'LineWidth',1,'DisplayName','Te','Color',.5*ones(1,3)) ;

% ax=gca;
% ax.YAxis.Label.FontSize=12;
% ax.YAxis.Label.FontWeight='normal';
% ax.YTick=[1:3];

%Add flat model:
% aux1=Y-(Y/U)*U;
% aux1=sqrt(mean(aux1.^2));%/sqrt(meanVar);
% aux1=conv(aux1,ones(1,binw)/binw,'valid'); %Smoothing
% plot(aux1,'LineWidth',1,'DisplayName','Flat','Color','k') ;


%Add data reproduce (05/03/2022)
% if lower==1
%     yhat=model{1}.Out;
% else

    
% end
legend('Location','NorthEastOutside','AutoUpdate','off')
yl=ax.YAxis.Limits;
pp=patch([40 350 350 40],[0 0 max(aux1) max(aux1)],.7*ones(1,3),'FaceAlpha',.5,'EdgeColor','none'); %PATR - PATS 
uistack(pp,'bottom')
ax.YAxis.Limits=yl;
axis tight
yticks('auto')
grid on
% ylim([0 1])
% set(gca,'YScale','log')

%% R^2 plot 

subplot(Nx,5,51:54)
hold on
yhat= model{1}.Out; %C * xhat' ; %yhat = C
% aux1 = 1 - sum((Y2- yhat).^2)./sum((Y2- mean(Y2)).^2);
aux1 = 1 - sum((dd).^2)./sum((Y- mean(Y)).^2);
aux2 = my_Rsquared_coeff(Y,yhat,1);
aux1=conv(aux1,ones(1,binw)/binw,'valid'); %Smoothing
plot(aux1,'LineWidth',2,'DisplayName',[num2str(size(C,2)), 'states'],'Color','r') ;


% aux2= 1 - sum((Y2- NNMF).^2)./sum((Y2- mean(Y2)).^2);
% aux2=conv(aux2,ones(1,binw)/binw,'valid'); %Smoothing
% plot(aux2,'LineWidth',2,'DisplayName','NNMF','Color','k') ;

% aux3= 1 - sum((Y2- Y_2).^2)./sum((Y2- mean(Y2)).^2);
% aux3=conv(aux3,ones(1,binw)/binw,'valid'); %Smoothing
% plot(aux3,'-.','LineWidth',2,'DisplayName','2 States','Color','b') ;
ylabel({'R^{2}'})
grid on
ax.YAxis.Label.FontSize=12;
legend('Location','NorthEastOutside','AutoUpdate','off')
pp=patch([40 350 350  40],[0 0 1 1],.7*ones(1,3),'FaceAlpha',.5,'EdgeColor','none'); %PATR - PATS 
uistack(pp,'bottom')
% axis tight
% ylim([0.8 1.05])
yticks('auto')

% subplot(Nx,Ny,3+9*Ny)
% hold on 
% 
% % aux1 = 1 - sum((Y2- yhat).^2)./sum((Y2- mean(Y2)).^2);
% aux1=conv(cross,ones(1,binw)/binw,'valid'); %Smoothing
% plot(aux1,'LineWidth',2,'DisplayName','3 states','Color','r') ;
% 
% aux2=conv(cross2,ones(1,binw)/binw,'valid'); %Smoothing
% plot(aux2,'--','LineWidth',2,'DisplayName','PCA','Color','k') ;
% 
% aux3=conv(cross3,ones(1,binw)/binw,'valid'); %Smoothing
% plot(aux3,'-.','LineWidth',2,'DisplayName','2 states','Color','b') ;
% ylabel({'Pearson'; 'Correlation (r)'})
% 
% grid on
% ax.YAxis.Label.FontSize=12;
% legend('Location','NorthEastOutside','AutoUpdate','off')
% yline(nanmean(aux1(1:40)))
% pp=patch([40 480 480 40],[0 0 1 1],.7*ones(1,3),'FaceAlpha',.5,'EdgeColor','none'); %PATR - PATS 
% uistack(pp,'bottom')
% axis tight
% yticks('auto')

%% %TURN ON FOR AFTEREEFECTS 
% subplot(Nx,5,56:59)
% % [pp,cc,aa]=pca(dd');%,'Centered','off');
% [pp,cc,aa]=pca((dd(:,481:681)'),'Centered','off');
% % [coeff,score,latent,tsquared,explained,mu]=pca(Yasym,'Centered','off');
% %%Input has to be row observation  and columns variables 
% %%pp -  Corresponding matrix of eigenvectors 
% %%cc - data projected in the principal component 
% %%aa - vector of eigent values 
% hold on
% aux1=conv(cc(:,1)',ones(1,binw)/binw,'valid');
% plot(aux1,'LineWidth',1) ;
% title('First PC of residual, mov. avg.')
% grid on
% hold on 
% subplot(Nx,5,[55 60])
% imagesc(reshape(pp(:,1),12, length(ytl))')
% colormap(flipud(map))
% caxis([-1 1])
% set(gca,'XTick',[],'YTick',yt,'YTickLabel',ytl,'FontSize',fs)
%%

%subplot(Nx,Ny,3+9*Ny)
%hold on
%aux1=conv2(Y,[-.5,1,-.5],'valid'); %Y(k)-.5*(y(k+1)+y(k-1));
%aux1=sqrt(sum(aux1.^2))/sqrt(1.5);
%aux1=aux1./(8.23+sqrt(sum(Y(:,2:end-1).^2))); %
%aux1=conv(aux1,ones(1,binw)/binw,'valid'); %Smoothing
%plot(aux1,'LineWidth',1) ;
%title('Instantaneous normalized std of data')
%grid on
%set(gca,'YScale','log')

end
%% Save fig
fName='OpenSans';
txt=findobj(gcf,'Type','Text');
set(txt,'FontName',fName);
ax=findobj(gcf,'Type','Axes');
set(ax,'FontName',fName);
for i=1:length(ax)
    ax(i).Title.FontWeight='normal';
end
