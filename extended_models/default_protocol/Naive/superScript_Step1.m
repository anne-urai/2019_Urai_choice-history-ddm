function [] = superScript_Step1(sub,iter)


%This file is compiled into a standalone application, using Matlab mcc,
%The compiled fiels are later on called by a shell script to facilitate
%paralle fitting in a cluster

% This function can also be used as is (in Matlab). It takes two string
% arguments: sub, is the subject's id and iter is an integer used to set
% the random seed. 

if ischar(sub), sb = str2double(sub); end
if ischar(iter), iters   = str2double(iter); end
load motionEnergyData_AnkeMEG.mat % see these *mat files in the "data_processing" folder
load descriptives.mat

QQ(isnan(QQ))=0;

 rand('seed',sb+iters);
 randn('seed',sb+iters);
    %%
coh=unique(data.behavior.coherence);

opts= optimset('Display','off','TolX',1e-4,'TolFun',1e-4,'MaxIter',150,'MaxFunEvals',1000); %options of simplex
drfts=((data.motionenergy_normalized(:,16:end))); % take motion energy samples after discarding the first 16
indx=find(data.behavior.subj_idx==sb & data.behavior.RT>0 & ~isnan(data.behavior.response) & ~isnan(data.behavior.prevresp)); %select relevant data
drfts=drfts(indx,:);
drfts=mean(drfts,2)/100; % take the mean of the motion energy samples and divided it by 100

matinput=[drfts data.behavior.coherence(indx) data.behavior.stimulus(indx) double(data.behavior.stimulus(indx)==data.behavior.prevresp(indx))];
mattinput=repmat(matinput,10,1); % define mattinput that contains informaiton relevat for the fitting. 
                                 % repeat each experimental trial 10x ti
                                 % improve robustness to simulation noise
                                 % see cost_fit for description of the
                                 % mattinput matrix

for tt=1:20 %generate 20 random starting points and keep the best as starting point to simplex
    stp=[rand(1,5)];
    out=cost_fit(sb,QQ,QN,mattinput,stp); %cost_fit is a function that calculates the -sum log likelihood given a set of parameters
    pars(tt,:)=stp;
    outp(tt)=out;
  
end
[k,l]=min(outp);
stp=pars(l,:);
[xbest,fx,exitflag,output]=fminsearchbnd(@(pars) cost_fit(sb,QQ,QN,mattinput,pars),[stp],[0 0.01 0 0 0],[1 1 1 1 1],opts);
% run simplex (fmnisearchbnd) and keep the best fitting parameter set for
% this particular subject and random seed (determined by iter).
out1=xbest;
out2=fx;

clear data;
 save(sprintf('/home/ktsetsos/code1/DDM/fits/out%s_%s.mat', sub, iter)); % saving path