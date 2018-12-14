function [] = AAM_Opt2(sub,iter)

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

opts= optimset('Display','off','TolX',1e-4,'TolFun',1e-4,'MaxIter',100,'MaxFunEvals',1000);

drfts=data.motionenergy_normalized(:,16:end);
drt=[];
for m=1:30;
    drt=[drt drfts(:,m) drfts(:,m) drfts(:,m) drfts(:,m) drfts(:,m)];
end
drfts=drt;
indx=find(data.behavior.RT>0 & data.behavior.subj_idx==sb  & ~isnan(data.behavior.response) & ~isnan(data.behavior.prevresp));

drfts=drfts(indx,:);
drfts=(drfts)/100;


out=NaN;
matinput=[0*data.behavior.coherence(indx) data.behavior.coherence(indx) data.behavior.stimulus(indx) double(data.behavior.stimulus(indx)==data.behavior.prevresp(indx))];
mattinput=repmat(matinput,15,1);
drfts=repmat(drfts,15,1);



load(sprintf('fits2/out%s.mat', sub))

tic;
mm=size(parms,1)
for tt=1:mm
tt
    stp=[parms(tt,:)];
    Rout(tt)=cost_fit(sb,QQ,QN,drfts,mattinput,stp);
end
Rout
toc;
[ifx,ixx]=sort(Rout);

tic;
for mm=1:5
stp=parms(ixx(mm),:);
[xbest,fx,exitflag,output]=fminsearchbnd(@(pars) cost_fit(sb,QQ,QN,drfts,mattinput,pars),[stp],[0 0.01 0 0 0],[1 1 1 1 1],opts);
Rout1(mm,:)=xbest;
Rout2(mm)=fx;
end
toc
clear data drfts drt mattinput matinput;
save(sprintf('/home/ktsetsos/code1/DDM/fits3/out%s_%s.mat', sub, iter));