function [] = superScript_Step1(sub,iter)

%%SAME AS THE CORRESPONDING SCRIPT IN THE NAIVE MODEL
%%DIFFERENCE 1:  IN LINE 31 WHERE THE mattinput matrix has an extra column
%%to store the previous response

%% DIFFERENCE 2: Instead of using 20 random starting points the simplex starting parameter set
%% takes the best fitting parameters from the Naive model and adds '0' for the extra bias parameter

if ischar(sub), sb = str2double(sub); end
if ischar(iter), iters   = str2double(iter); end
load motionEnergyData_AnkeMEG.mat
load descriptives.mat
load(sprintf('/home/ktsetsos/code1/DDM/fits/out%s_%s.mat',sub,iter)); 
stp2=xbest;
QQ(isnan(QQ))=0;

 rand('seed',sb+iters);
 randn('seed',sb+iters);
    %%
coh=unique(data.behavior.coherence);

opts= optimset('Display','off','TolX',1e-4,'TolFun',1e-4,'MaxIter',150,'MaxFunEvals',1000);
drfts=((data.motionenergy_normalized(:,16:end)));

indx=find(data.behavior.subj_idx==sb & data.behavior.RT>0 & ~isnan(data.behavior.response) & ~isnan(data.behavior.prevresp));

drfts=drfts(indx,:);
drfts=mean(drfts,2)/100;


out=NaN;
matinput=[drfts data.behavior.coherence(indx) data.behavior.stimulus(indx) double(data.behavior.stimulus(indx)==data.behavior.prevresp(indx)) data.behavior.prevresp(indx)];
mattinput=repmat(matinput,10,1);


    stp=[stp2 0];%starting parameter set is the best fit from the Naive model and 0 for the bias parameter
    out=cost_fit(sb,QQ,QN,mattinput,stp);


[xbest,fx,exitflag,output]=fminsearchbnd(@(pars) cost_fit(sb,QQ,QN,mattinput,pars),[stp],[0 0.01 0 0 0 -1],[1 1 1 1 1 1],opts);
out1=xbest;
out2=fx;

clear data;
 save(sprintf('/home/ktsetsos/code1/DDM_DC/fits/out%s_%s.mat', sub, iter));