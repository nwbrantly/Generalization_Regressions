% Script to separate the conditions 
%% Load params file 
subID='PATR06'
load([subID, 'params.mat'])

%%

%Adaptation 
load([subID, 'params.mat'])
adaptData = AddingConditions(adaptData, 'Adaptation', 'Adaptation2', true, 'Tied to Split');
save([subID, 'params.mat'],'adaptData')
changeCondName(subID,{'Adaptation', 'Adaptation2'},{'TM base','Adaptation'})
%%
% Tied to Pos Short 
load([subID, 'params.mat'])
adaptData = AddingConditions(adaptData, 'Pos short', 'Pos short 2', true, 'Tied to Split');
save([subID, 'params.mat'],'adaptData')
changeCondName(subID,{'Pos short', 'Pos short 2'},{'Tied Pos', 'Pos short'})



%Pos Short  to tied
load([subID, 'params.mat'])
adaptData = AddingConditions(adaptData, 'Pos short', 'Tied post Pos', false, 'Tied to Split');
save([subID, 'params.mat'],'adaptData')

%Tied to Neg Short 
load([subID, 'params.mat'])
adaptData = AddingConditions(adaptData, 'Neg short', 'Neg short 2', true, 'Tied to Split');
% adaptData = AddingConditions(adaptData, 'Tied Neg', 'Neg short', true, 'Tied to Split');
save([subID, 'params.mat'],'adaptData')
changeCondName(subID,{'Neg short', 'Neg short 2'},{'Tied Neg', 'Neg short'})

%Neg Short  to tied
load([subID, 'params.mat'])
adaptData = AddingConditions(adaptData, 'Neg short', 'Tied post Neg', false, 'Tied to Split');
save([subID, 'params.mat'],'adaptData')
%%
%Tied to Pos ramp  
% load('PATR04params.mat')
% adaptData = AddingConditions(adaptData, 'Pos short ramp', 'Pos short ramp 2', true, 'Tied to Split');
% save('PATR04params','adaptData')
% changeCondName('PATR04',{'Pos short ramp', 'Pos short ramp 2'},{'Tied Pos ramp', 'Pos short ramp'})

%Neg Short  to tied
load([subID, 'params.mat'])
adaptData = AddingConditions(adaptData, 'Pos short ramp', 'Tied post ramp', false, 'Tied to Split');
save([subID, 'params.mat'],'adaptData')

%% 
adaptData = AddingConditions(adaptData, 'Multiple pos short splits', 'Split Pos 1', true, 'Tied to Split');


% adaptData = AddingConditions(adaptData, 'Pos short ramp', 'Tied post ramp', false, 'Split to tied');
%%

for t= 1:30 
    if rem(t, 2) == 0
        adaptData = AddingConditions(adaptData, ['Tied Pos ' num2str(t-1)], ['Split Pos ' num2str(t)], true, ['Tied to Split', num2str(t)]);
        
    else
        if t==1
            adaptData = AddingConditions(adaptData, ['Split Pos ' num2str(t)], ['Tied Pos ' num2str(t)], false, ['Split to Tied', num2str(t)]);
        else
            adaptData = AddingConditions(adaptData, ['Split Pos ' num2str(t-1)], ['Tied Pos ' num2str(t)], false, ['Split to Tied', num2str(t)]);
        end
    end
end
save([subID, 'params.mat'],'adaptData')