figure


for steps=1:680
    
    if steps<40
        figure(1)
        surface(Casym(:,1),Casym(:,2),Yasym(steps,:))%,'filled', 'MarkerFaceColor','r')
        
        % title('Data')
    elseif steps>40 && steps<480
        figure(2)
        scatter3(Casym(:,1),Casym(:,2),Yasym(steps,:),'filled', 'MarkerFaceColor','g')
        ylabel('C_{contextual}')
        xlabel('C_{reactive}')
        zlabel('Yasym')
        % title('Data')
    elseif steps>480
        figure(3)
        scatter3(Casym(:,1),Casym(:,2),Yasym(steps,:),'filled', 'MarkerFaceColor','b')
        ylabel('C_{contextual}')
        xlabel('C_{reactive}')
        zlabel('Yasym')
        % title('Data')
    end
    hold on
    %     pause(.5)
end
% ylabel('C_{contextual}')
% xlabel('C_{reactive}')
% zlabel('Yasym')
% title('Data')
figure(1)
title('Baseline')
ylabel('C_{contextual}')
xlabel('C_{reactive}')
zlabel('Yasym')
set(gcf,'color','w')
figure(2)
title('Adaptation')
ylabel('C_{contextual}')
xlabel('C_{reactive}')
zlabel('Yasym')
set(gcf,'color','w')
figure(3)
title('Post-adapt')
ylabel('C_{contextual}')
xlabel('C_{reactive}')
zlabel('Yasym')
set(gcf,'color','w')

% postadapt=Yasym(41:480,:);
% [pp,cc,aa]=pca(postadapt,'Centered','off');


%%
figure


for steps=1:680
    scatter3(pp(:,1),pp(:,2),Yasym(steps,:),'filled')
    hold on 
%     pause(.5)
end
ylabel('PC1')
xlabel('PC2')
zlabel('Yasym')
title('PCA')
set(gcf,'color','w')