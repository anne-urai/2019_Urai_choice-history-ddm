function [] = superScript_Step2(sub,iter)

%%SAME AS THE CORRESPONDING SCRIPT IN THE NAIVE MODEL
%%DIFFERENCE:  IN LINE 31 WHERE THE mattinput matrix has an extra column
%%to store the previous response


if ischar(sub), sb = str2double(sub); end
if ischar(iter), iters   = str2double(iter); end
load motionEnergyData_AnkeMEG.mat
load descriptives.mat

QQ(isnan(QQ))=0;

 rand('seed',sb+iters);
 randn('seed',sb+iters);
    %%
coh=unique(data.behavior.coherence);

opts= optimset('Display','off','TolX',1e-4,'TolFun',1e-4,'MaxIter',100,'MaxFunEvals',1000);
drfts=((data.motionenergy_normalized(:,16:end)));

indx=find(data.behavior.subj_idx==sb & data.behavior.RT>0 & ~isnan(data.behavior.response) & ~isnan(data.behavior.prevresp));

drfts=drfts(indx,:);
drfts=mean(drfts,2)/100;




matinput=[drfts data.behavior.coherence(indx) data.behavior.stimulus(indx) double(data.behavior.stimulus(indx)==data.behavior.prevresp(indx)) data.behavior.prevresp(indx)];
mattinput=repmat(matinput,20,1);



load(sprintf('fits2/out%s.mat', sub))

tic;
mm=size(parms,1);
for tt=1:mm
tt
    stp=[parms(tt,:)];
    Rout(tt)=cost_fit(sb,QQ,QN,mattinput,stp);
end
Rout
toc;
[ifx,ixx]=sort(Rout);

tic;
for mm=1:5
stp=parms(ixx(mm),:);
[xbest,fx,exitflag,output]=fminsearchbnd(@(pars) cost_fit(sb,QQ,QN,mattinput,pars),[stp],[0 0 0.01 0 0 -1],[1 1 1 1 1 1],opts);
Rout1(mm,:)=xbest;
Rout2(mm)=fx;
end
toc
clear data;
save(sprintf('/home/ktsetsos/code1/DDM_DC/fits3/out%s_%s.mat', sub, iter));