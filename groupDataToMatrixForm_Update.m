function [Y,Yasym,Ycom,U,Ubreaks,Ysum,Yindv,labels]=groupDataToMatrixForm_Update(subjIdx,fName,sqrtFlag)
%% Load real data:
EMGdata=h5read(fName,'/EMGdata');
SLA=h5read(fName,'/SLA');
speedDiff=h5read(fName,'/speedDiff');
breaks=h5read(fName,'/breaks');
labels=hdf5read(fName,'/labels');

%%
U=speedDiff;
Ubreaks=breaks;

%% Some pre-proc
% if nargin<1
%     subjIdx=2:16; %Excluding C01 only
% end
%%
muscPhaseIdx=1:size(EMGdata,2);
Y=EMGdata(:,muscPhaseIdx,subjIdx);
% Yindv=nanmedian(Y,3);
Yindv=Y;

if nargin>2 && sqrtFlag
    Y=sqrt(Y);
end
Yasym=Y-fftshift(Y,2);
Ycom=Y-Yasym;
Yasym=Yasym(:,1:size(Yasym,2)/2,:);
Ysum=Y+fftshift(Y,2);
Ysum=Ysum(:,1:size(Ysum,2)/2,:);
Ycom=Ycom(:,1:size(Ycom,2)/2,:);
Y=nanmedian(Y,3); %Median across subjs
end
