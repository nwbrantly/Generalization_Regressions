%% Color definition
% Red, white and blue color display
%this patriotic color combination was decided to make the plot more
%accesible for colorblind people 

ex2=[0.2314    0.2980    0.7529];
ex1=[0.7255    0.0863    0.1608];
% ex1=[1,0,0];
% ex2=[0,0,1];
cc=[0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350      0.0780    0.1840];

mid=ones(1,3);
N=100;
% gamma=1.5; %gamma > 1 expands the white (mid) part of the map, 'hiding' low values. Gamma<1 does the opposite
gamma=1;
map=[flipud(mid+ (ex1-mid).*([1:N]'/N).^gamma); mid; (mid+ (ex2-mid).*([1:N]'/N).^gamma)];
