function [LL]= cost_fit(sub,QQ,QN,mattinput,pars)

coh1=[0  0.0900    0.2700]*100;
coh2=[0.03 0.0900  0.2700]*100;

for re=0:1
for c=1:length(coh1)
    
    trls=mattinput((mattinput(:,2)==coh1(c) | mattinput(:,2)==coh2(c)) & mattinput(:,4)==re,1);
    cor=mattinput((mattinput(:,2)==coh1(c)  | mattinput(:,2)==coh2(c)) & mattinput(:,4)==re,3);
    stim=mattinput((mattinput(:,2)==coh1(c)  | mattinput(:,2)==coh2(c)) & mattinput(:,4)==re,5); %this is the only difference relative
    %to the naive cost_fit function. Stim encodes the category of the
    %previous choice: +1 right, -1 left.

    [RT,ch]= diffusion_custom2(trls,cor,stim,[pars]);
   
    
    acc(c)=mean(ch);%*0.99+0.005;
    P(re+1,c,1,:)=(histcounts(RT(ch==0),[-inf squeeze(QQ(sub,re+1,c,1,:))' inf])./sum(ch==0))*(1-acc(c));
    P(re+1,c,2,:)=(histcounts(RT(ch==1),[-inf squeeze(QQ(sub,re+1,c,2,:))' inf])./sum(ch==1))*(acc(c));
    
   
end
end
P=P+eps;
NL=squeeze(QN(sub,:,:,:,:)).*log(P);
LL=-sum(sum(sum(sum(NL))));
