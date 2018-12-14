function [LL]= cost_fit(sub,QQ,QN,mattinput,pars)

% This function takes as inputs the subject number, the empirical 
% RT quantiles and QN (empirical trial count within each quantile)
% mattinput (see below) and the parameter vector
% It returns the sum of the negative log likelihood

%mattinput is trials X INFO matrix
% the INFO columns are: 1) mean motion energy, 2) coherence level, 3) stimulus
% identity, 4) stimulus consistency with previous choice (0 inconsistent, 1
% consistent).

coh1=[0  0.0900    0.2700  ]*100; % the 3 relevant coherence levels
coh2=[0.03 0.0900    0.2700  ]*100;

for re=0:1 % loop through STIM inconsistent/ consistent with previous choice
for c=1:length(coh1) % loop through coherence levels
    

    trls=mattinput((mattinput(:,2)==coh1(c) | mattinput(:,2)==coh2(c)) & mattinput(:,4)==re,1); % take relevant mean motion energy
    cor=mattinput((mattinput(:,2)==coh1(c)  | mattinput(:,2)==coh2(c)) & mattinput(:,4)==re,3); % take relevant stimulus identity

    [RT,ch]= diffusion_custom2(trls,cor,[pars]); %here the function implemented the model is called
                                                 % returns RT's and choices
                                                 % for all trials.
    acc(c)=mean(ch); % mean accuracy
    P(re+1,c,1,:)=(histcounts(RT(ch==0),[-inf squeeze(QQ(sub,re+1,c,1,:))' inf])./sum(ch==0))*(1-acc(c));
    P(re+1,c,2,:)=(histcounts(RT(ch==1),[-inf squeeze(QQ(sub,re+1,c,2,:))' inf])./sum(ch==1))*(acc(c));
    
    % We calculate the probability (P) predicted by the model for each empirical quantile (QQ) bin
   
end
end
P=P+eps; % avoid Nan's and Inf's in the logarithmic transformation
NL=squeeze(QN(sub,:,:,:,:)).*log(P); % Maximum likelihood over quantiles
LL=-sum(sum(sum(sum(NL))));