%3D betas representation

groupID={'ATS','ATR'};

load(['GroupRegression01_',groupID{1},'.mat'])

figure
ellipsoid(trans1.Coefficients.Estimate(1),trans1.Coefficients.Estimate(2),trans1.Coefficients.Estimate(3),trans1.Coefficients.SE(1),trans1.Coefficients.SE(2),trans1.Coefficients.SE(3))
hold on
if groupID=='ATR'
    plot3(0.8,0.2,0,'o','Color','#0072BD','MarkerFaceColor','#0072BD','MarkerSize',10)
    axis([.5 1 0 0.5  -.2 0.2])
else
    plot3(0.5,0,0.7,'o','Color','#0072BD','MarkerFaceColor','#0072BD','MarkerSize',10)
end
xlabel('Adapt')
ylabel('Noadapt')
zlabel('ContexSwitch')
grid on
% 
legend('$\hat{Post}_1$','$Post_1$','Interpreter','latex')
title('Post_1')
set(gcf,'color','w');

figure
hold on
ellipsoid(trans2.Coefficients.Estimate(1),trans2.Coefficients.Estimate(2),trans2.Coefficients.Estimate(3),trans2.Coefficients.SE(1),trans2.Coefficients.SE(2),trans2.Coefficients.SE(3))
if groupID=='ATR'
    plot3(0.2,0.2,1,'o','Color','#EDB120','MarkerFaceColor','#EDB120','MarkerSize',10)
    axis([0 .5 -.5 .5  .5 1])
else
    plot3(0.5,0.2,1,'o','Color','#EDB120','MarkerFaceColor','#EDB120','MarkerSize',10)
end
xlabel('Adapt')
ylabel('Noadapt')
zlabel('ContexSwitch')
grid on
% 
legend('$\hat{Post}_2$','$Post_2$','Interpreter','latex')
title('Post_2')
set(gcf,'color','w');
%% Transition together

figure
ellipsoid(trans1.Coefficients.Estimate(1),trans1.Coefficients.Estimate(2),trans1.Coefficients.Estimate(3),trans1.Coefficients.SE(1),trans1.Coefficients.SE(2),trans1.Coefficients.SE(3))
hold on
if groupID=='ATR'
    plot3(0.8,0.2,0,'o','Color','#0072BD','MarkerFaceColor','#0072BD','MarkerSize',10)
else
    plot3(0.5,0,0.7,'o','Color','#0072BD','MarkerFaceColor','#0072BD','MarkerSize',10)
end

xlabel('Adapt')
ylabel('Noadapt')
zlabel('ContexSwitch')

hold on
ellipsoid(trans2.Coefficients.Estimate(1),trans2.Coefficients.Estimate(2),trans2.Coefficients.Estimate(3),trans2.Coefficients.SE(1),trans2.Coefficients.SE(2),trans2.Coefficients.SE(3))
if groupID=='ATR'
    plot3(0.2,0.2,1,'o','Color','#EDB120','MarkerFaceColor','#EDB120','MarkerSize',10)
else
    plot3(0.5,0.2,1,'o','Color','#EDB120','MarkerFaceColor','#EDB120','MarkerSize',10)
end
xlabel('Adapt')
ylabel('Noadapt')
zlabel('ContexSwitch')
grid on
axis([0 1 -.5 .5  -.5 1])
set(gcf,'color','w');
legend('$\hat{Post}_1$','$Post_1$','$\hat{Post}_2$','$Post_2$','Interpreter','latex')
