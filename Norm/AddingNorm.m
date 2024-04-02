function [out]=AddingNorm(groupData,subID,newLabelPrefix,groupID)
%Adding Norm to groupAdaptationData

% Load data and computes norms for the entire time courses
%This code will find Euclidean norms for the entire time courses
%Created by DMM0 5/2022

%Modified 4/2024 DMMO

% 1) Computing Stride by stride norm
% 2) Compute bias removed stride by stride norm
% 3) Computing norm per muscle


% TODO: The code only gets the indiviudal muscle norm value we needs to implement: 
%1) Indivual muscle asymmetry 
%2) Removing the bias per muscle 

%%  1) Computing Stride by stride norm

for idx = 1:numel(subID)
    data=[];
    temp=[];
    aux1=[];
    
    subjIdx = find(contains(groupData.ID, subID{idx}));

    
    
    
    if ~isempty(subjIdx)
        
        Subj = groupData.adaptData{subjIdx};

        
        for i = 1:numel(newLabelPrefix)
            DataIdx=find(cellfun(@(x) ~isempty(x),regexp(Subj.data.labels,['^' newLabelPrefix{i} '[ ]?\d+$'])));
            
            data=[data Subj.data.Data(:,DataIdx)];
            data(isnan(data))=0;

        end
        
        data(isnan(data))=0;
        dataAsym=data-fftshift(data,2);
        dataAsym=dataAsym(:,1:size(dataAsym,2)/2,:);
        temp(:,1)=vecnorm(data');
        temp(:,2)=vecnorm(dataAsym');
        aux1=find(temp(:,1)>50);
        temp(aux1,:)=nan;
        groupData.adaptData{idx}.data=groupData.adaptData{idx}.data.appendData(temp,{'NormEMG','NormEMGasym'},...
            {'Norm of all the muscles','Norm asym of all the muscles'});
    end
    
    
    
end

%% 2) Compute bias removed stride by stride norm

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

for idx = 1:numel(subID)
    data=[];
    temp=[];
    data3=[];
    data3asym=[];

    subjIdx = find(contains(groupData.ID, subID{idx}));

    
    
    if ~isempty(subjIdx)
        
        Subj = groupData.adaptData{subjIdx};

        
        for i = 1:numel(newLabelPrefix)
            DataIdx=find(cellfun(@(x) ~isempty(x),regexp(Subj.data.labels,['^' newLabelPrefix{i} '[ ]?\d+$'])));            
                        
            data=[data Subj.data.Data(:,DataIdx)];
            data(isnan(data))=0;

        end    
        trial=find(contains(Subj.data.labels, {'trial'}));
        tt=unique(Subj.data.Data(:,trial));
        for t=1:length(tt)
            zz=tt(t);
            aux2=[];
            aux3=[];
            if find(contains(Subj.data.trialTypes(zz), {'OG'} ))
                
                Idx = find(Subj.data.Data(:,trial)==zz);
                aux2=data(Idx,:)';
                data2= aux2-OGref(:,subjIdx);
                
                aux3=aux2-fftshift(aux2,1);
                aux3=aux3(1:size(aux3,1)/2,:,:);
                
                data2asym=aux3-OGrefasym(:,subjIdx);
                
                
                
                
                
            else
                
                Idx = find(Subj.data.Data(:,trial)==zz);
                aux2=data(Idx,:)';
                data2= aux2-TRref(:,subjIdx);
                
                aux3=aux2-fftshift(aux2,1);
                aux3=aux3(1:size(aux3,1)/2,:,:);
                
                data2asym=aux3-TRrefasym(:,subjIdx);
                
                
            end
   
            data3=[data3 data2];
            data3asym=[data3asym data2asym];
            
        end
        data3(isnan(data3))=0;
        data3asym(isnan(data3asym))=0;
        temp(:,1)=vecnorm(data3);
        temp(:,2)=vecnorm(data3asym);
        aux1=find(temp(:,1)>50);
        temp(aux1,:)=nan;
        aux1=groupData.adaptData{idx}.data.Data;
        groupData.adaptData{idx}.data=groupData.adaptData{idx}.data.appendData(temp,{'UnBiasNormEMG','UnBiasNormEMGasym'},{'Context specifci unbais Norm of all the muscles','Context specifci unbais Norm asym of all the muscles'});
    end
    
    
    
end


%% 3) Computing norm per muscle
label=strcat(newLabelPrefix,'Norm');
desc=strcat(strcat(strcat(label,' muscle during the full gait cycle')));

for idx = 1:numel(subID)
    data=[];
    temp=[];
%     aux1=[];
    
    subjIdx = find(contains(groupData.ID, subID{idx}));

    if ~isempty(subjIdx)
        
        Subj = groupData.adaptData{subjIdx};
        
        
        for i = 1:numel(newLabelPrefix)
            
            DataIdx=find(cellfun(@(x) ~isempty(x),regexp(Subj.data.labels,['^' newLabelPrefix{i} '[ ]?\d+$'])));
           
            data=[Subj.data.Data(:,DataIdx)];
            data(isnan(data))=0;
            
            
            data(isnan(data))=0;
%             dataAsym=data-fftshift(data,2);
%             dataAsym=dataAsym(:,1:size(dataAsym,2)/2,:);
            temp(:,i)=vecnorm(data');
%             temp2(:,i)=vecnorm(dataAsym');
%             aux1=find(temp(:,1)>50);
%             temp(aux1,:)=nan;
        end
        
        groupData.adaptData{idx}.data=groupData.adaptData{idx}.data.appendData(temp,label,...
            desc);
    end
  
end

out= groupData;


end