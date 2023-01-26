%% Script to get the changes in EMG during Early and Late epochs (Adaptation and Post-adapt)
clear; close all; clc;
%% General variables

eE=1; %Number of strides early effect
eL=5; %Number of strides late
nstrides=5;% do not change

baseEp=defineEpochs({'Base'},{'TM base'}',[-40],[eE],[eL],'nanmean');
%% Adaptaton 
groupID ='BAT';
[normalizedGroupData, newLabelPrefix]=creatingGroupdataWnormalizedEMG(groupID);
%% Removing bad muscles 
%This script make sure that we always remove the same muscle for the
%different analysis 
normalizedGroupData= RemovingBadMuscleToSubj(normalizedGroupData);
% [normalizedGroupData]=RemoveBadMuscles(normalizedGroupData,{'BATS02','BATS04','BATS06','BATS09','BATS12'},{{'fSOLs','sSOLs','fVMs','sVMs','fVLs','sVLs','sRFs','fRFs'},{'fBFs','sBFs'},{'fRFs','sRFs'},{'fRFs','sRFs'},{'sRFs','fRFs','fVLs','sVLs'}});

%% Plotting 
ep=defineEpochs({'eA','lA'},{'Adaptation','Adaptation'},[ nstrides -40],[eE eE ],[eL eL],'nanmean');
% ep=defineEpochs({'EarlyAdapt'},{'Adaptation'},[nstrides],[eE],[eL],'nanmean');
% ep=defineEpochs({'EarlyPost'},{'Post 1'},[nstrides],[eE],[eL],'nanmean');
% baseEp=defineEpochs({'Base'},{'TM base'}',[-40],[eE],[eL],'nanmean');
fi=figure('Name',[groupID,' ','EMG structure']);
ph=tight_subplot(1,length(ep),[.03 .005],.04,.04);
summFlag='nanmedian';
fdr=.04;
savefig=0;
EMGchanges_Adaptation(normalizedGroupData,newLabelPrefix,ep,fi,ph,baseEp,fdr,summFlag,groupID,savefig) 


%% Post-adaptation 
%% Overgound 
% clear all
groupID ='BATS';
[normalizedGroupData, newLabelPrefix]=creatingGroupdataWnormalizedEMG(groupID);
%% Removing bad muscles 
%This script make sure that we always remove the same muscle for the
%different analysis 
normalizedGroupData= RemovingBadMuscleToSubj(normalizedGroupData);
%% Plotting 
ep=defineEpochs({'eP'},{'Post 1'},[nstrides],[eE],[eL],'nanmean');
% ep=defineEpochs({'PosShort','Ramp','eA','lA','eP','exp'},{'Pos short','Pos short ramp','Adaptation','Adaptation','Post 1','Split Pos 18'},...
%     [-10 -10 nstrides -40 nstrides -40],[eE eE eE eE eE eE],[eL eL eL eL eL eL],'nanmean');
% baseEp=defineEpochs({'Base'},{'TM base'}',[-40],[eE],[eL],'nanmean');
fi=figure('Name',[groupID,' ','EMG structure']);
ph=tight_subplot(1,length(ep),[.03 .005],.04,.04);
summFlag='nanmedian';
fdr=.05;
savefig=0;
EMGchanges_Adaptation(normalizedGroupData,newLabelPrefix,ep,fi,ph,baseEp,fdr,summFlag,groupID,savefig) 

%% Treadmill
% clear all
groupID ='BATR';
[normalizedGroupData, newLabelPrefix]=creatingGroupdataWnormalizedEMG(groupID);
%% Removing bad muscles 
%This script make sure that we always remove the same muscle for the
%different analysis 
normalizedGroupData= RemovingBadMuscleToSubj(normalizedGroupData);
%% Plotting
% ep=defineEpochs({'eA','lA','eP'},{'Adaptation','Adaptation','Post 1'},[ nstrides -40 5],[eE eE eE],[eL eL eL],'nanmean');
ep=defineEpochs({'eP'},{'Post 1'},[nstrides],[eE],[eL],'nanmean');
% ep=defineEpochs({'PosShort','Ramp','eA','lA','eP','exp'},{'Pos short','Pos short ramp','Adaptation','Adaptation','Post 1','Split Pos 18'},...
%     [-10 -10 nstrides -40 nstrides -40],[eE eE eE eE eE eE],[eL eL eL eL eL eL],'nanmean');
% baseEp=defineEpochs({'Base'},{'TM base'}',[-40],[eE],[eL],'nanmean');
fi=figure('Name',[groupID,' ','EMG structure']);
ph=tight_subplot(1,length(ep),[.03 .005],.04,.04);
summFlag='nanmedian';
fdr=.05;
savefig=0;
EMGchanges_Adaptation(normalizedGroupData,newLabelPrefix,ep,fi,ph,baseEp,fdr,summFlag,groupID,savefig) 

