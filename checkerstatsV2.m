function [p,h,alphaAdj]=checkerstatsV2(data1,data2,wsflag,parflag,alpha,posthoc,mindif)
%FUNCTION made by Digna deKam - use on deKam et al 2020
%function that computes stats on checkerboards
%Input:
%-data1 and 2: data matrices to compare
%-wsflag: 1 for within subjects comparison and 0 for between subjects
% comparison
%-parflag:0 for non-parametric, 1 for parametric
%-alpha level
%-posthoc test: default is 'benhoch'


if nargin<6
    posthoc='benhoch';
end

if nargin<7
    mindif=0;
end

if ndims(data1)==4
    data1=squeeze(data1);
    data2=squeeze(data2);
end
p=nan(size(data1,2),size(data1,1));



%perform actual test
for i=1:size(p,1)
    for k=1:size(p,2)
        if wsflag==0%perform between subject stats   
            if parflag==0
                p(i,k)=ranksum(squeeze(data1(k,i,:)),squeeze(data2(k,i,:)));         
            elseif parflag==1      
                [dummy,p(i,k)]=ttest2(squeeze(data1(k,i,:)),squeeze(data2(k,i,:)));
            end
        elseif wsflag==1%perform within subject stats    
            if parflag==0         
               % p(i,k)=signtest(squeeze(data1(k,i,:)),0); 
               p(i,k)=signrank(squeeze(data1(k,i,:)));%mindif has to be zero for two-tailed test  
                
               %pR(i,k)=signrank(squeeze(data1(k,i,:)),mindif,'tail','right','method','exact');         
               %pL(i,k)=signrank(squeeze(data1(k,i,:)),-mindif,'tail','left','method','exact');     
               %p=2*min(pR,pL);      
               %disp('performing ws stats')
                   
               
            elseif parflag==1        
                [dummy,p(i,k)]=ttest(squeeze(data1(k,i,:)),mindif);
                %[dummy,pR(i,k)]=ttest(squeeze(data1(k,i,:)),mindif,'tail','right');
                %[dummy,pL(i,k)]=ttest(squeeze(data1(k,i,:)),-mindif,'tail','left');
                %p=2*min(pR,pL); 
            end
        end
    end      
        
end


%perform posthoc corrections
if strcmp(posthoc,'none')
    h=zeros(size(p));
    h(find(p<alpha))=1;
    alphaAdj=alpha;
elseif strcmp(posthoc,'benhoch')
    ptemp=p(:);
    [htemp,alphaAdj,i1] = BenjaminiHochbergNew(ptemp,alpha,'true');
    h=reshape(htemp,size(p));    
elseif strcmp(posthoc,'bonferroni')
    h=zeros(size(p));
    alphaAdj=alpha/(length(p(:)));
    h(find(p<alphaAdj))=1;
end
    
    
    
end