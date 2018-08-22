
% ================================================= % 
% code from https://www.mathworks.com/matlabcentral/fileexchange/31661-fit-glm-with-quadratic-penalty?s_tid=prof_contriblnk
% ================================================= % 

clear all; close all;
addpath('/Users/urai/Documents/code/fitglmqp130');
path = '~/Data/psychophysicalKernels';
load(sprintf('%s/%s', path, 'motionEnergyData_AnkeMEG.mat'));

% grab data
r = double(data.behavior.response > 0);
X = [data.motionenergy_normalized * 0.01 ones(length(r), 1)];

% prior motion energy traces, separately for the final baseline
qf = blkdiag(qfsmooth1D(size(X, 2)-2), 0.01, 0.01);

% set some default options
opts.family   = 'binomlogit'; 
opts.Display  = 'final';

% 1. fit on all observers
results_fixedeffects  = cvglmfitqp(r,X,qf,getcvfolds(length(r), 50), opts);
figure; subplot(221); plot(data.timeaxis(1:end-1), results_fixedeffects.w(1:end-2)); 

% ================================================= % 
% NOW DO ACROSS OBSERVERS
% ================================================= % 

% opts.lambda0    = results_fixedeffects.finallambda; % what exactly does lambda0 do?

% now fit to each observer
results         = splitapply(@(r,X) cvglmfitqp(r,X, ...
    blkdiag(qfsmooth1D(size(X, 2)-2), 0.01, 0.01), ...
    getcvfolds(length(r), 10), opts), ...
    r, X, findgroups(data.behavior.transitionprob));

results_nocrossval = splitapply(@(r,X) glmfitqp(r,X, ...
    results_fixedeffects.finallambda * blkdiag(qfsmooth1D(size(X, 2)-2), 0.01, 0.01), opts), ...
    r, X, findgroups(data.behavior.subj_idx));

 disp(results);

% overview of the results
subplot(222); bar([results.finallambda]);
hline(results_fixedeffects.finallambda);
ylabel('Lambda'); xlabel('Subjects');

subplot(223); hold on;
weights = cat(2, results.w)';
weights = weights - weights(:, 1);

lambdas = [results.finallambda];
weights(lambdas > results_fixedeffects.finallambda * 10, :) = NaN;
for i = 1:size(weights, 1)
    plot(data.timeaxis(1:end-1), weights(i, 1:end-2), 'color', [1 0 lambdas(i)/max(lambdas)])
end
hold on;
boundedline(data.timeaxis(1:end-1), nanmean(weights(:, 1:end-2)), ...
    nanstd(weights(:, 1:end-2)) ./ sqrt(32), 'cmap', [0 0 0]);
%ylim([0.17 0.22]);
title('Individual lambda');

% also do the fit while enforcing a single lambda
subplot(224); hold on;
weights = cat(2, results_nocrossval.w)';
for i = 1:size(weights, 1)
    plot(data.timeaxis(1:end-1), weights(i, 1:end-2), 'color', [1 0 lambdas(i)/max(lambdas)])
end
hold on;
boundedline(data.timeaxis(1:end-1), mean(weights(:, 1:end-2)), ...
    std(weights(:, 1:end-2)) ./ sqrt(32), 'cmap', [0 0 0]);
%ylim([0.17 0.22]); 
title('Same lambda');

%% also make a plot that's the same as normal kernels
biasedkernels = results_fixedeffects.w(1:end-2)';
data.timeaxis = data.timeaxis(1:end-1);

% then average over coherence levels within each subject
close all; 
subplot(441);
hold on;
colors(1, :) = [0.3 0.3 0.3]; c = 1;
plot(data.timeaxis, (biasedkernels), 'color', colors(c, :), 'linewidth', 0.5);
plot(data.timeaxis(13:end), (biasedkernels(:, 13:end)), 'color', 'k', 'linewidth', 3);
axis tight;

ylabel({'Excess motion'; 'energy fluctuations (%)'});
xlabel('Time from stimulus onset (s)');
axis tight; xlim([0 0.75]); set(gca, 'xtick', 0:0.25:0.75);
box off; offsetAxes;
set(gca, 'xcolor', 'k', 'ycolor', 'k');
tightfig;
print(gcf, '-dpdf', '~/Data/serialHDDM/psychophysicalKernels_constrainedLogRes.pdf');

