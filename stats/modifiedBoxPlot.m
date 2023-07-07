function stats=modifiedBoxPlot(x,y,labels)
%Functon created by MAG to mimic the barplots from Wang and ﻿Srinivasan 2014 

%% First compute the statistics of each dataset

means = nanmean(y);
SDs = nanstd(y);
medians = nanmedian(y);
stats=[];
%% Plot the rectangles which will be located in mean +/- SD

% figure();
hold on;

for i=1:length(x)
    
       %Plot the box around the mean, since the box will be a SD from the mean
    rectangle('position',[x(i)-0.25 means(i)-SDs(i) 0.5 2*SDs(i)], 'FaceColor', 'w');
    xtemp = [x(i)-0.25:0.01:x(i)+0.25];
   
    %Plot the 95% confidence interval of the data
%     disp(['mean = ', num2str(means(i)),' CI ='])
     P = prctile(y(:,i),[2.5 97.5],"all");
   errorbar(x(i),means(i),abs(means(i)-P(1)), abs(means(i)-P(2)), 'Color', 'k');
%     errorbar(x(i),means(i),abs(P(1)), abs(P(2)), 'Color', 'k');



 %Plot mean and median
 plot(xtemp, means(i)*ones(length(xtemp),1), 'LineWidth', 1.5, 'Color', 'b');
 plot(xtemp, medians(i)*ones(length(xtemp),1), 'LineWidth', 1.5, 'Color', 'r', 'LineStyle', '--');

    %     %Plot the values that are outside the 95% interval
    idx = find(y(:,i)<P(1) | y(:,i)>P(2));
%    plot(x(i), y(idx,i),'Marker','+', 'Color', 'r');

stats.mean(i)=means(i);
stats.CI(:,i)=P;
end

newX = [x(1)-1 : 0.01 : x(end)+1];
plot(newX, zeros(length(newX),1), 'k');
hold off;
xlabel(labels)

end