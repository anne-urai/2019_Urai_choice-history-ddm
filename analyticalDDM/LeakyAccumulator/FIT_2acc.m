% function [best_ps,best_L,BIC] = FIT_2acc(subj,constraints)
%
% Fits the two-accumulator LCA to single-subject adaptive timescale data,
% using intial rough grid search and subsequent SIMPLEX optimization.
%
% Input: subj = string variable representing participant ID
%        constraints = string of parameter constraints (see main script)
%
% Outputs: 'best_ps' = matrix of best-fitting parameters per condition
%          'best_L' = likelihood value maximized by optimization procedure
%          'BIC' = BIC of best fit given n free parameters

function [best_ps,best_L,BIC,optim_ps,nLL] = FIT_2acc(subj,constraints)

% Adding required paths
addpath(genpath('D:\Experiments\Adaptive_timescale\Analysis'))
addpath(genpath('D:\Experiments\Adaptive_timescale\Analysis\FMINSEARCHBND'))

% Define load/save directories
loadpath = 'D:\Experiments\Adaptive_timescale\Task\My_task\Data\Behaviour\';
savepath = 'D:\Experiments\Adaptive_timescale\Analysis\Modelling\Basic2accumulator\Fits\';
if ~exist([savepath,subj,filesep,constraints],'dir')
    mkdir(savepath,[subj,filesep,constraints])
end
savepath = [savepath,subj,filesep,constraints,filesep];

% Choosing parameter constraints (switch b/w 'fi' for fixed, 'fr' for free and 'no' for none (leak only))
set_p.beta = constraints(1:2);            % brightness exponent
set_p.sigmaA = constraints(4:5);          % accumulator noise
set_p.theta = constraints(7:8);           % decision bound
set_p.leak = constraints(10:11);          % accumulator leak

disp('')
disp(['beta=',set_p.beta,', sigma=',set_p.sigmaA,', theta=',set_p.theta,', leak=',set_p.leak])

range_p.beta = [0.05 1];    % parameter ranges
range_p.sigmaA = [0 0.5];
range_p.theta = [0.1 inf];
range_p.leak = [0 0.5];

grid_p.beta = [0.2:0.2:1];    % parameter values over which to conduct initial rough grid search
grid_p.sigmaA = [0:0.03:0.21 0.26 0.31];
grid_p.theta = [0.4:0.06:0.94];  % NB - for theta, these are percentiles of max noise-only accumulator values, not absolute parameter values (which will vary with leak-noise pairs)
grid_p.leak = [0:0.02:0.18 0.2 0.25 0.3];

nsims.hit = 10000;  % number of model simulations per cell of experimental design
nsims.fa = 10000;

nOpt = 10;    % number of parameter sets, after grid search, to be passed on to Simplex routine
nIter = 500;  % number of Simplex iterations per parameter set

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Pulling behavioural data %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
allsubj = {'EIv' 'MHo' 'JRu' 'JHa' 'MSo' 'LMe' 'LHe' 'SJe' 'FSt' 'ARa' 'FRa' 'ROr' 'HKr' 'JMo' 'ZFe' 'NLe' 'TSc' 'SEn' 'HSa' 'DLe' 'TSu' 'RWi' 'URa' 'FMe' 'JFi' 'LBr' 'JSt' 'TSt' 'NKu' 'JNe'};
sessions = [4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4];
n_blocks = [12 11 12 12 12 12 12 11 12 12 12 12 11 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12;...
            12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12;...
            12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12;...
            12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12];
durorder = [1 2 1 2 1 1 1 1 2 2 2 1 2 2 1 1 1 1 2 2 2 1 2 2 1 1 2 2 1 2;...
            2 2 1 1 2 1 2 1 1 2 1 2 1 1 2 1 2 2 2 1 2 2 1 1 2 2 2 1 1 2;...
            1 1 2 1 2 2 1 2 1 1 1 2 2 2 2 2 1 2 1 1 1 2 2 2 1 1 1 2 2 1;...
            2 1 2 2 1 2 2 2 2 1 2 1 1 1 1 2 2 1 1 2 1 1 1 1 2 2 1 1 2 1];  % 1 = S-blocks first, 2 = L-blocks first
n_trials = 62;
Durs = [0.15 0.45 0.9];

discard_blocks = 1;   % number of blocks @ start of each run to discard (might take a while to adapt to statistics of new S/L condition)
discard_sessions = [0]; % whole sessions to discard (to potentially mitigate practice effects)

subji = find(strcmp(allsubj,subj));
disp(sprintf('Pulling behavioural data for subject %d of %d: %s',subji,length(allsubj),subj))
Behav_all_L=[]; Behav_all_S=[];
for s = 1:sessions(subji);
    if sum(ismember(discard_sessions,s))==0
        
        if (strcmp(subj,'JNe') && s==1) || (strcmp(subj,'JMo') && s==1) % dumping 1st n blocks for a couple of subjects due to changing difficulty in 1st session
            startb = 4;
        else startb = 1;
        end
        
        for b = startb:n_blocks(s,subji);
            load([loadpath,subj,filesep,'S',num2str(s),filesep,subj,'_',num2str(s),'_',num2str(b),'.mat']) % load behaviour
            Behav = Behav(1:n_trials,:);  % pull only relevant trials (early mistake led to error here for first participant or two such that extra ones for these subejcts should be ignored)
            
            if (b<=6 && b>discard_blocks) || (b>6+discard_blocks)
                if durorder(s,subji) == 1
                    if b<=6
                        Behav_all_S(end+1:end+n_trials,1:10) = Behav;
                    else Behav_all_L(end+1:end+n_trials,1:10) = Behav;
                    end
                else
                    if b<=6
                        Behav_all_L(end+1:end+n_trials,1:10) = Behav;
                    else Behav_all_S(end+1:end+n_trials,1:10) = Behav;
                    end
                end
            end
        end
    end
end

Strs = sort(unique(Behav_all_L(find(Behav_all_L(:,4)>0),4)));  % pulling signal strengths
if length(Strs)>2  % In case there's more than two signal strengths, pulling two most common ones
    for i=1:length(Strs);
        nStrs(i) = length(find(Behav_all_L(:,4)==Strs(i)));
    end
    [~,inds] = sort(nStrs);
    Strs = Strs(inds(end-1:end));
end

% Calculating accuracy for each condition/trial-type
data.fa(1,1:2) = calc_input(Behav_all_S,51,0,[]);
data.fa(2,1:2) = calc_input(Behav_all_L,51,0,[]);

for s = 1:length(Strs);
    for d = 1:length(Durs);
        data.hit(1,s,d,1:2) = calc_input(Behav_all_S,1,Durs(d),Strs(s));
        data.hit(2,s,d,1:2) = calc_input(Behav_all_L,1,Durs(d),Strs(s));
    end
end

% Storing number of observations for BIC calculation
nObs = size(Behav_all_S,1)+size(Behav_all_L,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Running rough grid search %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist([savepath,'grid_results.mat'])  % if grid search has not been performed yet
    % Loop through all grid combinations of beta, sigmaA and leak to determine appropriate threshold values;
    % create matrix of full parameter sets covering entire grid; calculate likelihoods as we go
    gridSets=[]; gridL=[];
    disp(sprintf('Starting rough grid search'))
    if ~strcmp(set_p.leak,'no')   % if leak is included
        for b = 1:length(grid_p.beta);
            disp(sprintf('%d / %d of grid searched...',b-1,length(grid_p.beta)))
            for s = 1:length(grid_p.sigmaA);
                for l = 1:length(grid_p.leak);
                    colvals = (size(gridSets,1)+1):(size(gridSets,1)+length(grid_p.theta));  % getting positions into which to add new parameter sets
                    gridSets(colvals,1) = grid_p.beta(b);    % adding this round of betas
                    gridSets(colvals,2) = grid_p.sigmaA(s);  %       "      "       sigmas
                    gridSets(colvals,4) = grid_p.leak(l);    %       "      "       leaks
                    
                    y = sim2acc_no_bound([grid_p.beta(b) grid_p.sigmaA(s) grid_p.leak(l) 0.6],nsims.fa);   % simulating noise-only trials
                    gridSets(colvals,3) = quantile(max(max(y,[],3)),grid_p.theta);  % calculating appropriate threshold values based on maximum accumulator values given this beta-noise-leak combo
                    
                    for t = 1:length(colvals);
                        for n = 1:size(data.fa,1);  % looping through conditions to calculate FA likelihoods given current parameter set
                            gridL(colvals(t),n) = L2acc(data.fa(n,:),1-grid_p.theta(t));
                        end
                        
                        for str = 1:length(Strs);  % looping through signal strengths/durations and running simulations for this parameter set
                            for dur = 1:length(Durs);
                                ACC = sim2acc([gridSets(colvals(t),:) 0.6],Strs(str),Durs(dur),nsims.hit);
                                for n = 1:size(data.fa,1);  % looping through conditions to calculate HIT likelihoods given current parameter set
                                    gridL(colvals(t),n) = gridL(colvals(t),n)*L2acc(data.hit(n,str,dur,:),length(find(ACC==1))/length(ACC));
                                end
                            end
                        end
                    end
                end
            end
        end
    else   % if leak is fixed @ 0
        for b = 1:length(grid_p.beta);
            disp(sprintf('%d / %d of grid searched...',b-1,length(grid_p.beta)))
            for s = 1:length(grid_p.sigmaA);
                colvals = (size(gridSets,1)+1):(size(gridSets,1)+length(grid_p.theta));  % getting positions into which to add new parameter sets
                gridSets(colvals,1) = grid_p.beta(b);    % adding this round of betas
                gridSets(colvals,2) = grid_p.sigmaA(s);  %       "      "       sigmas
                
                y = sim2acc_no_bound([grid_p.beta(b) grid_p.sigmaA(s) 0 0.6],nsims.fa);   % simulating noise-only trials
                gridSets(colvals,3) = quantile(max(max(y,[],3)),grid_p.theta);  % calculating appropriate threshold values based on maximum accumulator values given this beta-noise-leak combo
                
                for t = 1:length(colvals);
                    for n = 1:size(data.fa,1);  % looping through conditions to calculate FA likelihoods given current parameter set
                        gridL(colvals(t),n) = L2acc(data.fa(n,:),1-grid_p.theta(t));
                    end
                    
                    for str = 1:length(Strs);  % looping through signal strengths/durations and running simulations for this parameter set
                        for dur = 1:length(Durs);
                            ACC = sim2acc([gridSets(colvals(t),:) 0 0.6],Strs(str),Durs(dur),nsims.hit);
                            for n = 1:size(data.fa,1);  % looping through conditions to calculate HIT likelihoods given current parameter set
                                gridL(colvals(t),n) = gridL(colvals(t),n)*L2acc(data.hit(n,str,dur,:),length(find(ACC==1))/length(ACC));
                            end
                        end
                    end
                end
            end
        end
    end
    disp('Grid search complete...')
    
    % Saving grid search results
    save([savepath,'grid_results.mat'],'data','gridSets','gridL','grid_p','set_p')
else
    load([savepath,'grid_results.mat'])
end

% Recover best n parameter sets from set of viable parameter combinations given constraints (***********CURRENTLY ONLY WORKS WHEN FITTING TWO CONDITIONS!!!!!!***********)
[max_ps,max_Ls] = get_grid_min(gridSets,gridL,grid_p,set_p,nOpt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Running finer parameter optimization %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Creating matrix of parameter constraints (n x p: n = number of conditions, p = number of parameters per condition; when element=0, this parameter is fixed @ zero)
par_mat = get_param_mat(set_p,size(data.fa,1));

% Constructing parameter starting point vectors
init_ps = [];
for s = 1:length(max_ps);
    current_set=[];
    for p = 1:size(max_ps{s},2);
        current_set(end+1:end+length(unique(par_mat(:,p)))) = max_ps{s}(1:unique(par_mat(:,p)),p);
    end
    init_ps(s,1:length(current_set)) = current_set;
end

% Constructing parameter range matrix
range_p_rep{1}=range_p.beta; range_p_rep{2}=range_p.sigmaA; range_p_rep{3}=range_p.theta; range_p_rep{4}=range_p.leak; 
par_range=[];
for p = 1:size(par_mat,2);
    if sum(par_mat(:,p))>0  % making sure this parameter is actually being used
        par_range(end+1:end+length(unique(par_mat(:,p))),1:2) = repmat(range_p_rep{p},length(unique(par_mat(:,p))),1);
    end
end

% Looping through all starting points and executing optimization
options = optimset('fminsearch');  % specifying Simplex routine settings
options.MaxIter = nIter;  % maximum number of iterations

disp(sprintf(''))
for s = 1:size(init_ps,1);
    disp(sprintf('Starting SIMPLEX optimization %d of %d...',s,size(init_ps,1)))
    [optim_ps(s,:),nLL(s,1)] = fminsearchbnd(@(pm) ML_2acc(pm,par_mat,data,Strs,Durs,nsims),init_ps(s,:),par_range(:,1)',par_range(:,2)',options);  % fitting
end

best_ps_v = optim_ps(find(nLL==min(nLL)),:);  % pull best-fitting parameter set
best_L = nLL(find(nLL==min(nLL)));

% Convert best-fitting params from vector into matrix with rows=conditions
best_ps = zeros(size(par_mat));
for p = 1:size(par_mat,2);  % counting number of free params per parameter type
    if sum(unique(par_mat(:,p)))>0
        free_ps(p) = length(unique(par_mat(:,p)));
    else free_ps(p) = 0;  % setting to zero if parameter is fixed @ zero
    end
end
for p = 1:length(free_ps);  % contructing parameter sets for each distinct condition
    if free_ps(p)==0  % if parameter is fixed @ zero
        best_ps(:,p) = 0;
    else
        for f = 1:free_ps(p);
            if p==1
                best_ps(find(par_mat(:,p)==f),p) = best_ps_v(f);
            else best_ps(find(par_mat(:,p)==f),p) = best_ps_v(f+sum(free_ps(1:p-1)));
            end
        end
    end
end

% Print results
disp(sprintf('\nBest-fitting parameters:'))
disp(sprintf('S-duration: b=%1.3f, sigma=%1.3f, theta=%1.3f, leak=%1.3f',best_ps(1,:)))
disp(sprintf('L-duration: b=%1.3f, sigma=%1.3f, theta=%1.3f, leak=%1.3f',best_ps(2,:)))
disp(sprintf('Negative log likelihood: %1.1f',best_L))

% Calculate BIC of best-fitting model
BIC = (2*best_L)+(length(best_ps_v)*log(nObs));

% Saving best fits after optimization
save([savepath,'best_fit.mat'],'data','set_p','best_ps','best_L','BIC','init_ps','optim_ps','nLL')