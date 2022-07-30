%speed profiles update (Exploration)
clear all;close all;clc
mid=[0.5 1];
small=mid;% [0.65 0.85];

large=mid;%[0.35 1.15];

split={small, mid, large};
speed=[0.75*ones(50,2)];

for t=1:3
y=randsample(3,3);

speed= [speed; split{y(1)}.*ones(50,2);0.75*ones(randsample(10:30,1),2);...
    split{y(2)}.*ones(50,2);0.75*ones(randsample(10:30,1),2);...
    split{y(3)}.*ones(50,2);0.75*ones(randsample(10:30,1),2)];


end

speed=[speed; mid.*ones(150,2)];
%%
% figure
% plot(speed(:,1))
% hold on 
% plot(speed(:,2))

velR=speed(:,2);
velL=speed(:,1);

figure
plot(velR)
hold on
plot(velL)


%% 
legend('R','L')

%% Shor plit 

mid=[0.5 1];
speed=[];
speed=[0.75*ones(50,2); mid.*ones(30,2);0.75*ones(25,2)];

velR=speed(:,2);
velL=speed(:,1);

figure
plot(velR)
hold on
plot(velL)
legend('R','L')

% Negative 
NewR=velL;
NewL=velR;

velR=NewR;
velL=NewL;

figure
plot(velR)
hold on
plot(velL)
legend('R','L')
%% Adaptation 
mid=[0.5 1];
speed=[];
speed=[0.75*ones(50,2); mid.*ones(450,2)];
velR=speed(:,2);
velL=speed(:,1);

figure
plot(velR)
hold on
plot(velL)
legend('R','L')

%%

velR=[0.75*ones(50,1); velR];
velL=[0.75*ones(50,1); velL];

figure
plot(velR)
hold on
plot(velL)
legend('R','L')


