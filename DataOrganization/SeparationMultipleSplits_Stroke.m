% Script to separate the conditions
% In the EMG genralization protocol we did not stop the treadmill between
% transitions. Thus, we are separating this epochs in a post-hoc manner.
% this is a script intended for indvdal params files 
%% Load params file 
subID='C3S09_S1'
load([subID, 'params.mat'])



%% Short perturbations 
% Tied to Pos Short 
load([subID, 'params.mat'])
adaptData = AddingConditions(adaptData, 'Short Pos', 'Short Pos 2', true, 'Tied to Split');
save([subID, 'params.mat'],'adaptData')
changeCondName(subID,{'Short Pos', 'Short Pos 2'},{'Tied Pos', 'Short Pos'})

% Pos Short  to tied
load([subID, 'params.mat'])
adaptData = AddingConditions(adaptData, 'Short Pos', 'Tied post Pos', false, 'Tied to Split');
save([subID, 'params.mat'],'adaptData')

% Tied to Neg Short 
load([subID, 'params.mat'])
adaptData = AddingConditions(adaptData, 'Short Neg', 'Short Neg 2', true, 'Tied to Split');
save([subID, 'params.mat'],'adaptData')
changeCondName(subID,{'Short Neg', 'Short Neg 2'},{'Tied Neg', 'Short Neg'})

% Neg Short  to tied
load([subID, 'params.mat'])
adaptData = AddingConditions(adaptData, 'Short Neg', 'Tied post Neg', false, 'Tied to Split');
save([subID, 'params.mat'],'adaptData')




