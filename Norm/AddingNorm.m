function [out]=AddingNorm(groupData,subID,newLabelPrefix,groupID)
%Adding Norm to groupAdaptationData

% Computes norms for the entire time courses
%This code will find Euclidean norm the participants muscle activity for the entire time courses
%Created by DMM0 5/2022

%Modified 4/2024 DMMO

% 1) Computing Stride by stride norm
% 2) Compute bias removed stride by stride norm
% 3) Computing norm per muscle


% TODO: The code only gets the indiviudal muscle norm value we needs to implement: 
%1) Indivual muscle asymmetry 
%2) Removing the bias per muscle 

%%  1) Computing Stride by stride norm

for idx = 1:numel(subID) %Loop across participant 
    m_data=[]; %temp variables 
    temp=[];
    aux1=[];
    
    subjIdx = find(contains(groupData.ID, subID{idx})); %This can be use if we want to skip participants

    if ~isempty(subjIdx)
        
        Subj = groupData.adaptData{subjIdx}; 

        
        for i = 1:numel(newLabelPrefix)
            DataIdx=find(cellfun(@(x) ~isempty(x),regexp(Subj.data.labels,['^' newLabelPrefix{i} '[ ]?\d+$']))); %Finding the columns of the data 
            
            m_data=[m_data Subj.data.Data(:,DataIdx)]; %concatenating the data 
            m_data(isnan(m_data))=0; %nan are made zero to computer the norm 

        end
        
%         data(isnan(data))=0;
        dataAsym=m_data-fftshift(m_data,2); % EMGasym
        dataAsym=dataAsym(:,1:size(dataAsym,2)/2,:); %Only keeping the fist half of the data 
        temp(:,1)=vecnorm(m_data'); %computing the euclidean norm 
        temp(:,2)=vecnorm(dataAsym');
%         aux1=find(temp(:,1)>50); %This was used for noisy data. CAUTION
%         temp(aux1,:)=nan;
        groupData.adaptData{idx}.data=groupData.adaptData{idx}.data.appendData(temp,{'NormEMG','NormEMGasym'},...
            {'Norm of all the muscles','Norm asym of all the muscles'});
    end
    
    
    
end

%% 2) Compute bias removed stride by stride norm

%  reference data to remove 
if contains(groupID,'NTS') ||  contains(groupID,'NTR') ||  contains(groupID,'CTR') || contains(groupID,'CTS') 
    ep=defineEpochVR_OG_UpdateV8('nanmean');
    refEpTR= defineReferenceEpoch('TRbase',ep);
    refEpOG= defineReferenceEpoch('OGbase',ep);
else
    ep=defineEpochs_regressionYA('nanmean');
    refEpTR= defineReferenceEpoch('TM base',ep);
    refEpOG= defineReferenceEpoch('OG base',ep);
end

padWithNaNFlag=true;

[OGref]=groupData.getPrefixedEpochData(newLabelPrefix,refEpOG,padWithNaNFlag); 
OGref=squeeze(OGref);
OGrefasym=OGref-fftshift(OGref,1);
OGrefasym=OGref(1:size(OGref,1)/2,:,:);

[TRref]=groupData.getPrefixedEpochData(newLabelPrefix,refEpTR,padWithNaNFlag); 
TRref=squeeze(TRref);
TRrefasym=TRref-fftshift(TRref,1);
TRrefasym=TRref(1:size(TRref,1)/2,:,:);


for idx = 1:numel(subID) %Loop across participant  
    
    m_data=[];%temp variables 
    temp=[];
    unbiasDataAll=[];
    unbiasDataasymAll=[];

    subjIdx = find(contains(groupData.ID, subID{idx}));  %This can be use if we want to skip participants

    
    
    if ~isempty(subjIdx)
        
        Subj = groupData.adaptData{subjIdx};

        
        for i = 1:numel(newLabelPrefix)
            DataIdx=find(cellfun(@(x) ~isempty(x),regexp(Subj.data.labels,['^' newLabelPrefix{i} '[ ]?\d+$'])));      %Finding the columns of the data       
                        
            m_data=[m_data Subj.data.Data(:,DataIdx)]; %concatenating the data 
            m_data(isnan(m_data))=0; %nan are made zero to computer the norm 

        end    
        
        trial=find(contains(Subj.data.labels, {'trial'})); %finding the label trials 
        tt=unique(Subj.data.Data(:,trial)); %Finding trials 
        for t=1:length(tt)
            
            trial_number=tt(t);
            aux2=[];
            aux3=[];
            
            if find(contains(Subj.data.trialTypes(trial_number), {'OG'} )) %If the data type is OG remove OG baseline
                
                Idx = find(Subj.data.Data(:,trial)==trial_number); %fiding the indices where the type is OG
                aux2=m_data(Idx,:)'; %grabing data OG data
                unbiasData= aux2-OGref(:,subjIdx); % Removing bias
                
                aux3=aux2-fftshift(aux2,1); %computing EMGaysm
                aux3=aux3(1:size(aux3,1)/2,:,:); %savoing the top part of the matrix
                
                unbiasDataasym=aux3-OGrefasym(:,subjIdx);v % Removing EMGasym bias
                
                
                
                
                
            elseif find(contains(Subj.data.trialTypes(trial_number), {'TM'} )) %If the data type is OG remove OG baseline
                
                Idx = find(Subj.data.Data(:,trial)==trial_number);
                aux2=m_data(Idx,:)';
                unbiasData= aux2-TRref(:,subjIdx);
                
                aux3=aux2-fftshift(aux2,1);
                aux3=aux3(1:size(aux3,1)/2,:,:);
                
                unbiasDataasym=aux3-TRrefasym(:,subjIdx);
                
            else
                
                warning('Update code to match your protocol and conditions' )
                
            end
   
            unbiasDataAll=[unbiasDataAll unbiasData];
            unbiasDataasymAll=[unbiasDataasymAll unbiasDataasym];
            
        end
        unbiasDataAll(isnan(unbiasDataAll))=0;
        unbiasDataasymAll(isnan(unbiasDataasymAll))=0;
        temp(:,1)=vecnorm(unbiasDataAll);
        temp(:,2)=vecnorm(unbiasDataasymAll);
%         aux1=find(temp(:,1)>50); %Clip the norm to less than 50 
%         temp(aux1,:)=nan;
        aux1=groupData.adaptData{idx}.data.Data;
        groupData.adaptData{idx}.data=groupData.adaptData{idx}.data.appendData(temp,{'UnBiasNormEMG','UnBiasNormEMGasym'},{'Context specifci unbais Norm of all the muscles','Context specifci unbais Norm asym of all the muscles'});
    end
    
    
    
end


%% 3) Computing norm per muscle
label=strcat(newLabelPrefix,'Norm'); %Creating label for the params file 
desc=strcat(strcat(strcat(label,' muscle during the full gait cycle'))); %Creating description for params file 

for idx = 1:numel(subID) %loop across participant 
    m_data=[];
    temp=[];
%     aux1=[];
    
    subjIdx = find(contains(groupData.ID, subID{idx}));  %This can be use if we want to skip participants

    if ~isempty(subjIdx)
        
        Subj = groupData.adaptData{subjIdx};
        
        
        for i = 1:numel(newLabelPrefix) %loop across muscle 
            
            DataIdx=find(cellfun(@(x) ~isempty(x),regexp(Subj.data.labels,['^' newLabelPrefix{i} '[ ]?\d+$'])));  %Finding the columns of the data       
            m_data=[Subj.data.Data(:,DataIdx)];
            m_data(isnan(m_data))=0;
            
            
            m_data(isnan(m_data))=0; %changing nan to zerot to compute norm 
            temp(:,i)=vecnorm(m_data'); %individual muscle norm 
%             aux1=find(temp(:,1)>50); %Caution clips norm to <50
%             temp(aux1,:)=nan;
        end
        
        groupData.adaptData{idx}.data=groupData.adaptData{idx}.data.appendData(temp,label,...
            desc);
    end
  
end

out= groupData;


end