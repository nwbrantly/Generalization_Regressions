% Grouping muscle by joint bar plots 
load('BATR_post-adaptation_29-November-2023.mat')

[muscle,subj]=find(R2>=0.5); % At the momment of cuttoff values in 50% of variance.
%If the model can reproduce 50% of the variance then it is not a godo model 

sank=[];
ship=[];
sknee=[];
fank=[];
fhip=[];
fknee=[];

for l=1:length(muscle)
    
            if muscle(l)<=3 %'sHIP'
                if isempty(ship)==1
                    ship(1,:)=[mean(reactive_trace(41:45,muscle(l),subj(l)),'omitnan') mean(contextual_trace(41:45,muscle(l),subj(l)),'omitnan') muscle(l)];

                else
                    ship(end+1,:)=[ mean(reactive_trace(41:45,muscle(l),subj(l)),'omitnan') mean(contextual_trace(41:45,muscle(l),subj(l)),'omitnan') muscle(l)];
                end
                
            elseif muscle(l)>=4 &&  muscle(l)<=9 %sTHIGH
                if isempty(sknee)==1
                    sknee(1,:)=[mean(reactive_trace(41:45,muscle(l),subj(l)),'omitnan') mean(contextual_trace(41:45,muscle(l),subj(l)),'omitnan') muscle(l)];
                else
                    sknee(end+1,:)=[mean(reactive_trace(41:45,muscle(l),subj(l)),'omitnan') mean(contextual_trace(41:45,muscle(l),subj(l)),'omitnan') muscle(l)];
                end

                
            elseif muscle(l)>=10 &&  muscle(l)<=14 %sSHANK
                
                if isempty(sank)==1
                    sank(1,:)=[mean(reactive_trace(41:45,muscle(l),subj(l)),'omitnan') mean(contextual_trace(41:45,muscle(l),subj(l)),'omitnan') muscle(l)];
                else
                    sank(end+1,:)=[mean(reactive_trace(41:45,muscle(l),subj(l)),'omitnan') mean(contextual_trace(41:45,muscle(l),subj(l)),'omitnan') muscle(l)];
                end
                
                
            elseif muscle(l)>=15 &&  muscle(l)<=17 %fHIP
             if isempty(fhip)==1
                    fhip(1,:)=[mean(reactive_trace(41:45,muscle(l),subj(l)),'omitnan') mean(contextual_trace(41:45,muscle(l),subj(l)),'omitnan') muscle(l)];
                else
                    fhip(end+1,:)=[ mean(reactive_trace(41:45,muscle(l),subj(l)),'omitnan') mean(contextual_trace(41:45,muscle(l),subj(l)),'omitnan') muscle(l)];
                end

                
            elseif muscle(l)>=18 &&  muscle(l)<=23 %fTHIGH
                if isempty(fknee)==1
                    fknee(1,:)=[mean(reactive_trace(41:45,muscle(l),subj(l)),'omitnan') mean(contextual_trace(41:45,muscle(l),subj(l)),'omitnan') muscle(l) ];
                else
                    fknee(end+1,:)=[mean(reactive_trace(41:45,muscle(l),subj(l)),'omitnan') mean(contextual_trace(41:45,muscle(l),subj(l)),'omitnan') muscle(l)];
                end

            elseif muscle(l)>=24 %fSHANK
                if isempty(fank)==1
                    fank(1,:)=[mean(reactive_trace(41:45,muscle(l),subj(l)),'omitnan') mean(contextual_trace(41:45,muscle(l),subj(l)),'omitnan') muscle(l)];
                else
                    fank(end+1,:)=[mean(reactive_trace(41:45,muscle(l),subj(l)),'omitnan') mean(contextual_trace(41:45,muscle(l),subj(l)),'omitnan') muscle(l)];
                end

            end
end



%% % bar plot by muscle group 
close all
colors=["#77AC30";"#7E2F8E";"#D95319"];

for c=1:2
figure(c)
hold on 
bar([1 4], [mean(ship(:,c)) mean(fhip(:,c))],.2,'FaceColor',colors(1))
errorbar([1 4],[mean(ship(:,c)) mean(fhip(:,c))],[std(ship(:,c))/sqrt(length(ship)) std(fhip(:,c))/sqrt(length(ship))] ,'.k')
plot(1.1,ship(:,c),'*k')
plot(4.1,fhip(:,c),'*k')

bar([2 5], [mean(sknee(:,c)) mean(fknee(:,c))],.2,'FaceColor',colors(2))
errorbar([2 5],[mean(sknee(:,c)) mean(fknee(:,c))],[std(sknee(:,c))/sqrt(length(sknee)) std(fknee(:,c))/sqrt(length(sknee))] ,'.k')
plot(2.1,sknee(:,c),'*k')
plot(5.1,fknee(:,c),'*k')

hold on 
bar([3 6], [mean(sank(:,c)) mean(fank(:,c))],.2,'FaceColor',colors(3))
errorbar([3 6],[mean(sank(:,c)) mean(fank(:,c))],[std(sank(:,c))/sqrt(length(sank)) std(fank(:,c))/sqrt(length(sank)) ] ,'.k')
plot(3.1,sank(:,c),'*k')
plot(6.1,fank(:,c),'*k')

xlim([0 7])
set(gca,'XTick',[1:6],'XTickLabel',{'sHIP','sKNEE','sANK','fHIP','fKNEE','fANK'},'FontSize',10)
xline(3.5)

if c==1
    ylabel('W_{reactive}')
else
    ylabel('W_{contextual}')

end

end 



