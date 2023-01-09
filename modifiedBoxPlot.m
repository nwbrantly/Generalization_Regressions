function modifiedBoxPlot(x,y,labels)

%% First compute the statistics of each dataset

means = nanmean(y,1);
SDs = nanstd(y,[],1);
medians = nanmedian(y,1);

%% Plot the rectangles which will be located in mean +/- SD

figure();
hold on;

for i=1:length(x)
    
       %Plot the box around the mean, since the box will be a SD from the mean
    rectangle('position',[x(i)-0.25 means(i)-SDs(i) 0.5 2*SDs(i)], 'FaceColor', 'w');
    xtemp = [x(i)-0.25:0.01:x(i)+0.25];
   
    %Plot the 95% confidence interval of the data
     P = prctile(y(:,i),[2.5 97.5],"all");
   errorbar(x(i),medians(i),abs(medians(i)-P(1)), abs(medians(i)-P(2)), 'Color', 'k');



 %Plot mean and median
 plot(xtemp, means(i)*ones(length(xtemp),1), 'LineWidth', 1.5, 'Color', 'b');
 plot(xtemp, medians(i)*ones(length(xtemp),1), 'LineWidth', 1.5, 'Color', 'r', 'LineStyle', '--');

    %     %Plot the values that are outside the 95% interval
    idx = find(y(:,i)<P(1) | y(:,i)>P(2));
   plot(x(i), y(idx,i),'Marker','+', 'Color', 'r');

end

newX = [x(1)-1 : 0.01 : x(end)+1];
plot(newX, zeros(length(newX),1), 'k');
hold off;
% xlabel(labels)

end