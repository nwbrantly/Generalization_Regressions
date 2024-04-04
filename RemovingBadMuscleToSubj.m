function groupData= RemovingBadMuscleToSubj(groupData)
%Code to make sure we remove the same muscle for all the analysis

% For now we need to hard code the participants ID and the muscle that we
% want to remove 

%INPUT:
%GroupData class

%OUTPUT:
%GroupData class with bad muscles replaced with NaN

%%
% Remmoved muscle removed the entire muscle for the data. This muscle are
% being removed due to bad normalization. 
%You can see the subjects and the correspondant muscle that is being remove

 % Treadmill group (Experiment 1)
[groupData]=RemoveBadMuscles(groupData,{'BATR14','BATR03'},{{'fVLs','sVLs','fVMs','sVMs'},{'fVLs','sVLs','fVMs','sVMs'}});

% Generalization group (Experiment 2)
 [groupData]=RemoveBadMuscles(groupData,{'BATS02','BATS04','BATS06','BATS09','BATS12'},...
     {{'fSOLs','sSOLs','fVMs','sVMs','fVLs','sVLs','sRFs','fRFs'},{'fBFs','sBFs'},{'fRFs','sRFs'},{'fRFs','sRFs'},{'sRFs','fRFs','fVLs','sVLs'}});
 
% Nimbus testing (R01)
 [groupData]=RemoveBadMuscles(groupData,{'NTS_03','NTS_07'},...
     {{'sLGs','sSEMTs','sLGs','sSEMTs'},{'fRFs','sRFs'},{'fRFs','sRFs'}});
 

 % Nimbus testing (R01)
 [groupData]=RemoveBadMuscles(groupData,{'CTS_03','CTS_05','CTS_04'},...
     {{'fHIPs'},{'sLGs'},{'sLGs','fTAs'}});
 
  % Very artificial groups testing (R01)
 [groupData]=RemoveBadMuscles(groupData,{'VATR03','VATS03','VATS04','VATS06','VATR04','VATR05','VATR06','VATS05'},...
     {{'sHIPs'},{'fPERs','fHIPs'},{'sRFs','sVMs'},{'fVMs','fPERs'},{'sRFs'},{'sVMs'},{'sHIPs'},{'fVMs','fVLs','fRFs'}});

 
end
 

