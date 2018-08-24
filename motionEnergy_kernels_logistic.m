
% ================================================= % 
% code from https://www.mathworks.com/matlabcentral/fileexchange/31661-fit-glm-with-quadratic-penalty?s_tid=prof_contriblnk
% ================================================= % 

clear all; close all;
addpath('/Users/urai/Documents/code/fitglmqp130');
path = '~/Data/psychophysicalKernels';
load(sprintf('%s/%s', path, 'motionEnergyData_AnkeMEG.mat'));

% normalize to range of -1, 1 to be in line with the prior
scaleFactor                  = 0.004; % determined through trial and error... this value ensures that the fitted weights match glmfit [0.004]
data.motionenergy_normalized = data.motionenergy_normalized * scaleFactor;
data.behavior.motionenergy   = mean(data.motionenergy_normalized, 2);
data.behavior.stimulus       = data.behavior.stimulus .* data.behavior.coherence .* scaleFactor;

% remove the initial filter rise time
timeStart                    = 5; 
data.motionenergy_normalized = data.motionenergy_normalized(:, timeStart:end);
data.timeaxis                = data.timeaxis(timeStart:end);

% grab data
r = double(data.behavior.response > 0);
X = [ones(length(r), 1) data.behavior.stimulus data.behavior.motionenergy];

% set some default options
opts.family   = 'binomlogit'; % assume logistic regression

% first do a normal glmfit to know the range of expected weights across subjects
psychfunc = splitapply(@(r,X) {glmfit(X, r, 'binomial', 'constant', 0)}, ...
    r, X, findgroups(data.behavior.subj_idx));
results.w_psychfunc = cat(2, psychfunc{:})';
results.w_psychfunc_fixed = glmfit(X, r, 'binomial', 'constant', 0);

% MAKE SURE THIS CODE CAN FIT A BASIC PSYCHOMETRIC FUNCTION FOR EACH SUBJECT!
results_fit = splitapply(@(r,X) evidenceglmfitqp(r,X, ...
    blkdiag( results.w_psychfunc_fixed(1), results.w_psychfunc_fixed(2), results.w_psychfunc_fixed(2)), opts), ...
    r, X, findgroups(data.behavior.subj_idx));
results.w_psychfunc_evidence = [results_fit.w]';

% do these two methods give the same result?
figure;
for i = 1:3,
    subplot(3,3,i); scatter(results.w_psychfunc(:, i), results.w_psychfunc_evidence(:, i));
    rl = refline(1); rl.Color = 'k'; lsline;
    title(sprintf('r = %.2f', corr(results.w_psychfunc(:, i), results.w_psychfunc_evidence(:, i))));
    xlabel('glmfit'); ylabel('evidence fit');
end

%% ================================================= % 
% NOW FIT ACROSS OBSERVERS' REAL RESPONSES
% ================================================= % 

% match the prior values to the known slope
qf = blkdiag(results.w_psychfunc_fixed(3) * qfsmooth1D(size(data.motionenergy_normalized, 2)), ...
     results.w_psychfunc_fixed(1));
X = [data.motionenergy_normalized ones(length(r), 1)];

% confirm that the slope and bias terms estimated are similar
results_fit_fixedeffects            = evidenceglmfitqp(r, X, qf, opts);
results.w_psychfunc_kernelfit_fixedeffects = results_fit_fixedeffects.w(end);
results.w_kernelfit_fixedeffects    = results_fit_fixedeffects.w(1:end-1);

% then, for each subject
results_fit     = splitapply(@(r,X) evidenceglmfitqp(r,X,qf,opts), ...
    r, X, findgroups(data.behavior.subj_idx));
results.w       = cat(2, results_fit.w)';

% compute the psychometric slope as the average of weights over time
results.w_psychfunc_kernelfit = [nanmean(results.w(:, 1:end-1), 2) results.w(:, end)];
results.w       = results.w(:, 1:end-1);

% plot the results
clf; subplot(221); plot(data.timeaxis, results.w_kernelfit_fixedeffects, 'k'); ylabel('Fixed effects');
hold on; yyaxis right; plot(data.timeaxis, nanmean(results.w), 'linewidth', 2); 
ylabel('Average of subjects');
axis tight; vline(data.timeaxis(13));

% these are the individual results
subplot(222); hold on;
plot(data.timeaxis, results.w);

hold on;
boundedline(data.timeaxis, mean(results.w), ...
    std(results.w) ./ sqrt(size(results.w, 1)), 'cmap', [0 0 0]);
%ylim([0.17 0.22]); 

subplot(223);
scatter(results.w_psychfunc_evidence(:, 3), results.w_psychfunc_kernelfit(:, 1));
xlabel('psychfunc slope');
ylabel('kernel slope');
title(sprintf('r = %.3f', corr(results.w_psychfunc_evidence(:, 3), results.w_psychfunc_kernelfit(:, 1))));
lsline;

% check that the intercepts do match as a sanity check
assert(corr(results.w_psychfunc(:, 1), results.w_psychfunc_kernelfit(:, 2)) > 0.9, ...
     'logistic intercepts dont match');
 

 %% ================================================= % 
% FIRST, CONFIRM THAT WE GET SENSIBLE WEIGHTS WHEN WE KNOW WHAT THE
% OBSERVER DOES
% ================================================= % 

% % SIMULATE A LINEARLY DECAYING WEIGHT AND AN OFFSET
% % FINDING: THE MODEL CAN RECOVER WEIGHTS BETWEEN 2-1 BEST!
% w_simul = results.w_psychfunc_fixed(2) .* ones(size(data.timeaxis))';
% w_simul = linspace(results.w_psychfunc_fixed(2)+2, results.w_psychfunc_fixed(2)-1, ...
%     numel(data.timeaxis))';
% 
% % simulate weights of the same order of magnitude, is this retrieved?
% r_simul = X*[w_simul; results.w_psychfunc_fixed(1)];
% r_simul = binornd(1,1./(1+exp(-r_simul)));
% 
% results_simulated = splitapply(@(r,X) evidenceglmfitqp(r,X, ...
%     blkdiag(qfsmooth1D(size(X, 2)-1), 0.01), opts), ...
%     r_simul, X, findgroups(data.behavior.prevresp));
% sim.w = cat(2, results_simulated.w)';
% sim.w = sim.w(:, 1:end-1);
% 
% results_simulated_fixedeffects = evidenceglmfitqp(r_simul, X, ...
%     blkdiag(qfsmooth1D(size(X, 2)-1), 0.01), opts);
% sim.w_fixedeffects = results_simulated_fixedeffects.w(1:end-1);
% 
% figure; subplot(221); 
% plot(w_simul, 'k', 'linewidth', 3); hold on;
% plot(sim.w_fixedeffects', 'k--', 'linewidth', 3); hold on;
% plot(mean(sim.w), 'k:', 'linewidth', 3);
% plot(sim.w');
