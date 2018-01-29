% Calculates total likelihood (ML) of data given a set of model parameters,
% as determined by the 2-accumulator LCA model
%
% Inputs:
%   params      = model parameters (to be adjusted by estimation routine...)
%   param_mat   = matrix of parameter constraints
%   data        = number of signal-present/hit (data.hit) and noise-only/false
%                 alarm (data.fa) trials in observed data
%   strs        = levels of signal strength
%   durs        = levels of signal duration
%   nsims       = number of model simulations for signal-present (nsims.hit) 
%                 and noise-only (nsims.fa) trial-types

function ML = ML_2acc(params,param_mat,data,strs,durs,nsims)

% Pull params for each condition: pm(n,1:5) = [beta, sigmaA, theta, leak, Ter]
pm = zeros(size(param_mat));
for p = 1:size(param_mat,2);  % counting number of free params per parameter type
    if sum(unique(param_mat(:,p)))>0
        free_ps(p) = length(unique(param_mat(:,p)));
    else free_ps(p) = 0;  % setting to zero if parameter is fixed @ zero
    end
end
for p = 1:length(free_ps);  % contructing parameter sets for each distinct condition
    if free_ps(p)==0  % if parameter is fixed @ zero
        pm(:,p) = 0;
    else
        for f = 1:free_ps(p);
            if p==1
                pm(find(param_mat(:,p)==f),p) = params(f);
            else pm(find(param_mat(:,p)==f),p) = params(f+sum(free_ps(1:p-1)));
            end
        end
    end
end
pm(:,end+1) = 0.6;  % fixing Ter = 0.6
            
% Simulate behaviour given parameters and calculate likelihood terms for each condition separately
L = zeros(1,size(pm,1));
for n = 1:size(pm,1);
    % Noise-only trials
    tic
    ACC = sim2acc(pm(n,:),0,0,nsims.fa);  % simulate behaviour
    toc
    L(n) = L2acc(data.fa(n,:),length(find(ACC==51))/length(ACC));  % calculate likelihood term
    % Signal-present trials
    for s = 1:length(strs);
        for d = 1:length(durs);
            tic
            ACC = sim2acc(pm(n,:),strs(s),durs(d),nsims.hit);  % simulate behaviour
            toc
            L(n) = L(n)*L2acc(data.hit(n,s,d,:),length(find(ACC==1))/length(ACC));  % increment likelihood term
        end
    end
end

% Take negative log of product of all likelihoods (i.e. quantity to be minimized)
ML = -log(prod(L));

