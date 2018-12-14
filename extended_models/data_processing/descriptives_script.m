%% The aim of this script is to summarise the raw data into quantiles on a
%% per participant basis.


clear all;
close all;
load('motionEnergyData_AnkeMEG.mat'); % +200 ms are added to all RT's to address fitting
                                      %issues with fast responses. This
                                      %structure is used in the default
                                      %protocol models.
%load('motionEnergyData_AnkeMEG_Dynamic.mat');% RT's are memasured to
                                       %stimulus onset. Therefore +750 ms
                                       %is added to all RT's. This
                                       %structure is used in the dynamic
                                       %protocol.

ss=unique(data.behavior.subj_idx);
cohs=unique(data.behavior.coherence);
cohs=cohs(1:end);
cohs1=[0 9 27];% 81% coherence trials are excluded
cohs2=[3 9 27];% 0% and 3% trials are pooled together due to indistinguishable
               % behaviour in these trials.

drfts=(mean(data.motionenergy_normalized(:,13:end),2)); % calculate the mean motion energy
                                                        % after discarding the 13 first
                                                        % samples (rising  time of the
                                                        % motion energy  filter)
                                                 
data.behavior.correct=double(sign(drfts)==sign(data.behavior.response)); 
rep=double(data.behavior.stimulus==data.behavior.prevresp);%1 if stimulus category is the same as previous response



for s=1:length(ss);

    for re=0:1 % loop across rep=0,1
    for c=1:length(cohs1) % loop across all coherence levels
        acc(re+1,ss(s),c)=nanmean(data.behavior.correct(data.behavior.subj_idx==ss(s) & rep==(re) & data.behavior.RT>0 & (data.behavior.coherence==cohs1(c) | data.behavior.coherence==cohs2(c))));
        % calculate accuracy for STIM consistent with previous response and
        % each coherence level
        for cor=0:1 % calculate incorrect/ correct quantiles
    indx=find(data.behavior.subj_idx==ss(s) &  rep==(re) &  data.behavior.RT>0 & (data.behavior.coherence==cohs1(c) | data.behavior.coherence==cohs2(c)) & data.behavior.correct==cor);
    RR=data.behavior.RT(indx);RR(isnan(RR))=[];
    QQ(ss(s),re+1,c,cor+1,:)=quantile(RR,[.1 .3 .5 .7 .9]);% QQ->the RT quantiles
    a=quantile(RR,[.1 .3 .5 .7 .9]);
    QN(ss(s),re+1,c,cor+1,1)=sum(RR<=a(1)); %QN--> number of trials on the corresponding quantile
    QN(ss(s),re+1,c,cor+1,2)=sum(RR>a(1) & RR<=a(2));
    QN(ss(s),re+1,c,cor+1,3)=sum(RR>a(2) & RR<=a(3));
    QN(ss(s),re+1,c,cor+1,4)=sum(RR>a(3) & RR<=a(4));
    QN(ss(s),re+1,c,cor+1,5)=sum(RR>a(4) & RR<=a(5));
    QN(ss(s),re+1,c,cor+1,6)=sum(RR>a(5));
    
        end
    end
    end
end

save descriptives acc QQ QN