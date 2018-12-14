function [] = superScript_Step1(sub,iter)

%% LARGELY SAME AS THE CORRESPONDING FUNCTION IN THE DEFAULT PROTOCOL
%% THE IMPORTANT DIFFERENCE IS THAT THE MOTION ENERGY IS NOT AVERAGED
%% BUT ITS WHOLE TIME-COURSE IS FED TO THE MODEL 


if ischar(sub), sb = str2double(sub); end
if ischar(iter), iters   = str2double(iter); end
load motionEnergyData_AnkeMEG.mat
load descriptives.mat

QQ(isnan(QQ))=0;

 rand('seed',sb+iters);
 randn('seed',sb+iters);
    %%
coh=unique(data.behavior.coherence);

opts= optimset('Display','off','TolX',1e-4,'TolFun',1e-4,'MaxIter',150,'MaxFunEvals',1000);
drfts=data.motionenergy_normalized(:,16:end);
drt=[];
for m=1:30;
    drt=[drt drfts(:,m) drfts(:,m) drfts(:,m) drfts(:,m) drfts(:,m)]; % interpolate motion energy since the dt of the model is different from the dt in the motion filtering
end
drfts=drt;
indx=find(data.behavior.RT>0 & data.behavior.subj_idx==sb  & ~isnan(data.behavior.response) & ~isnan(data.behavior.prevresp));

drfts=drfts(indx,:); %store the motion energy time-series on a per trial basis
drfts=(drfts)/100; 


out=NaN;
matinput=[0*data.behavior.coherence(indx) data.behavior.coherence(indx) data.behavior.stimulus(indx) double(data.behavior.stimulus(indx)==data.behavior.prevresp(indx))];
mattinput=repmat(matinput,10,1);
drfts=repmat(drfts,10,1);
for tt=1:20
    stp=[rand(1,5)];
    out=cost_fit(sb,QQ,QN,drfts,mattinput,stp);
    pars(tt,:)=stp;
    outp(tt)=out;
  
end
[k,l]=min(outp)
stp=pars(l,:);
[xbest,fx,exitflag,output]=fminsearchbnd(@(pars) cost_fit(sb,QQ,QN,drfts,mattinput,pars),[stp],[0 0.01 0 0 0],[1 1 1 1 1],opts);
out1=xbest;
out2=fx;

clear data;
clear drt drfts mattinput matinput;
 save(sprintf('/home/ktsetsos/code1/DDM/fits/out%s_%s.mat', sub, iter));