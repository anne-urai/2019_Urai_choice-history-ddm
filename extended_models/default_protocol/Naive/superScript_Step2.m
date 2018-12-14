function [] = superScript_Step2(sub,iter)

%This file is also compiled into a standalone application, using Matlab mcc,
%The compiled fiels are later on called by a shell script to facilitate
%paralle fitting in a cluster

% This function can also be used as is (in Matlab). It takes two string
% arguments: sub, is the subject's id and iter is an integer used to set
% the random seed. 


if ischar(sub), sb = str2double(sub); end
if ischar(iter), iters   = str2double(iter); end
load motionEnergyData_AnkeMEG.mat %see these *mat files in the "data_processing" folder
load descriptives.mat

QQ(isnan(QQ))=0;

 rand('seed',sb+iters);
 randn('seed',sb+iters);
    %%
coh=unique(data.behavior.coherence);

opts= optimset('Display','off','TolX',1e-4,'TolFun',1e-4,'MaxIter',100,'MaxFunEvals',1000);
drfts=((data.motionenergy_normalized(:,16:end))); % take motion energy samples after discarding the first 16

indx=find(data.behavior.subj_idx==sb & data.behavior.RT>0 & ~isnan(data.behavior.response) & ~isnan(data.behavior.prevresp));
% choose relevant data
drfts=drfts(indx,:);
%drfts=abs(drfts).^expn.*sign(drfts);
drfts=mean(drfts,2)/100;

matinput=[drfts data.behavior.coherence(indx) data.behavior.stimulus(indx) double(data.behavior.stimulus(indx)==data.behavior.prevresp(indx))];
mattinput=repmat(matinput,10,1); % define mattinput that contains informaiton relevat for the fitting. 
                                 % repeat each experimental trial 20x ti
                                 % improve robustness to simulation noise
                                 % see cost_fit for description of the
                                 % mattinput matrix



load(sprintf('fits2/out%s.mat', sub))

tic;
mm=size(parms,1)
for tt=1:mm % loop through the fits calculated in Step1 (50 per subject)
tt
    stp=[parms(tt,:)]; %store the parameters
    Rout(tt)=cost_fit(sb,QQ,QN,mattinput,stp); % recalculate the goodness of fit for each of the 50 iterations. 
                                               % This will be slightly
                                               % different relative to
                                               % Step_1 because the random
                                               % seed is different
                                          
end
toc;
[ifx,ixx]=sort(Rout); % sort the 50 fits by goodness of fit

tic;
for mm=1:5 % take the 5 best ones and use them as starting point to simplex
stp=parms(ixx(mm),:);
[xbest,fx,exitflag,output]=fminsearchbnd(@(pars) cost_fit(sb,QQ,QN,mattinput,pars),[stp],[0 0.01 0 0 0],[1 1 1 1 1],opts);
Rout1(mm,:)=xbest; % save the simplex output
Rout2(mm)=fx;
end
toc
clear data;
save(sprintf('/home/ktsetsos/code1/DDM/fits3/out%s_%s.mat', sub, iter)); %save to path 

%% At the end these 5 best fits are re-evaluated using a mattinput of 30x size. The best out of these 5 fits
%% is deemed the best fitting parameter set.