function motionEnergy_kernels()
addpath('~/Desktop/code/gramm');

%{
    Psychophysical kernels were computed by averaging contrast fluctuations as a
function of the sample number of the test stimulus. We treated the reference grating
as if it had been shown at the same time of the test stimulus for the same duration.
The psychophysical kernel at time point t is then given by:
E(t) = hE(t)S i + hE(t)N i
E(t)S is the contrast fluctuation of the selected option at time t and E(t)N is the c
ontrast fluctuation of the non selected stimulus at time t. Absolute stimulus contrast
was transformed into contrast fluctuation by subtracting average generative contrast
values (e.g. QUEST threshold) for the respective trial type (0.5 + threshold or 0.5 ?
threshold). Average generative contrast values were computed in blocks of 100 trials,
corresponding to natural blocks in the experiment after which participants took a
short break. The expectation is across trials.
%}

close all;
path = '~/Data/psychophysicalKernels';

% either all the data or only neutral
load(sprintf('%s/%s', path, 'motionEnergyData_AnkeMEG.mat'));

% only select neutral
% data.motionenergy               = data.motionenergy([data.behavior.transitionprob] == 0.5, :);
% data.motionenergy_normalized    = data.motionenergy_normalized([data.behavior.transitionprob] == 0.5, :);
% data.behavior                   = data.behavior([data.behavior.transitionprob] == 0.5, :);


% use the normalized motion energy, so that the units are % coherence 'up'
data.motionenergy = data.motionenergy_normalized;

selectTrialNum = 0; % select based on each individual's number of trials?
if ~selectTrialNum
    % only select trials without any objective evidence
    data.motionenergy = data.motionenergy([data.behavior.coherence] < 10, :);
    data.behavior     = data.behavior([data.behavior.coherence] < 10, :);
end

% what is the time-course of evidence that leads subjects to make their
% preferred vs non-preferred choice?
kernelFun       = @(x, y) nanmean(x(y, :)) - nanmean(x(~y, :));

% recode into choices that are biased or not
data.behavior.repeat = (data.behavior.response == data.behavior.prevresp);

% for each observers, compute their bias
[gr, sjs]   = findgroups(data.behavior.subj_idx);
sjrep       = splitapply(@nanmean, data.behavior.repeat, gr);
sjrep       = sjs(sjrep < 0.5);

% recode into biased and unbiased choices
data.behavior.biased = data.behavior.repeat;
altIdx      = ismember(data.behavior.subj_idx, sjrep);
data.behavior.biased(altIdx) = double(~(data.behavior.biased(altIdx))); % flip

% bin
coh                   = [data.behavior.stimulus .* data.behavior.coherence];
[gr, sj, coh]         = findgroups(data.behavior.subj_idx, coh);
biasedkernels         = splitapply(kernelFun, data.motionenergy, (data.behavior.response > 0), gr);

% average within each subject!
biasedkernels = splitapply(@nanmean, biasedkernels, findgroups(sj));

% =============================== %
% PSYCHOPHYSICAL KERNELS - without bias
% =============================== %

% then average over coherence levels within each subject
subplot(441);
hold on;
colors(1, :) = [0.3 0.3 0.3]; c = 1;
plot(data.timeaxis, nanmean(biasedkernels), 'color', colors(c, :), 'linewidth', 0.5);
b{c} = boundedline(data.timeaxis(13:end), nanmean(biasedkernels(:, 13:end)), ...
    nanstd(biasedkernels(:, 13:end)) ./ sqrt(length(unique(sj))), 'cmap', colors(c, :), 'alpha');
plot(data.timeaxis(13:end), nanmean(biasedkernels(:, 13:end)), 'color', colors(c, :), 'linewidth', 1);
axis tight;

% do statistics on the timecourse
[h, p, stat] = ttest_clustercorr(biasedkernels);
%[h, pval] = ttest(biasedkernels);
%[h, crit_p] = fdr_bh(pval, 0.05);

% remove significance during filter rise time
h(1:12) = 0;

ylims = get(gca, 'Ylim');
mask = double(h);
mask(mask==0) = nan;
mask = ((ylims(2)*0.1)+ylims(1))*mask; % plot a tiny bit above the lower ylim
plot(data.timeaxis, mask, '.', 'MarkerSize', 10, 'color', 'k');

ylabel({'Excess motion'; 'energy fluctuations (%)'});
xlabel('Time from stimulus onset (s)');
axis tight; xlim([0 0.75]); set(gca, 'xtick', 0:0.25:0.75);
box off; offsetAxes;
set(gca, 'xcolor', 'k', 'ycolor', 'k');
tightfig;
print(gcf, '-dpdf', '~/Data/serialHDDM/psychophysicalKernels.pdf');


%% COMPARE THE KERNELS WITH THE O-U EFFECTIVE LEAK PARAMETER
global mypath;
results     = readtable(sprintf('%s/summary/%s/allindividualresults_kostis.csv', mypath, 'Anke_MEG_Neutral'));
lambda      = results.ouK_vanilla_dv;
kernelDiff  = nanmean(biasedkernels(:, 13:29), 2) - nanmean(biasedkernels(:, 30:45), 2);
[rho, pval] = corr(kernelDiff, lambda, 'type', 'spearman');
fprintf('\n\nCorrelation between kernels and O-U lambda, Spearmans rho: %.3f, p = %.3f \n', rho, pval);

% % =============================== %
% % PSYCHOPHYSICAL KERNELS
% % =============================== %
% 
% % average those within each subject?
% close all;
% plot(data.timeaxis, kernelFun(data.motionenergy, (data.behavior.response > 0)));
% 
% % bin
% coh                   = [data.behavior.stimulus .* data.behavior.coherence];
% [gr, sj, coh, pref]   = findgroups(data.behavior.subj_idx, coh, data.behavior.biased);
% biasedkernels         = splitapply(kernelFun, data.motionenergy, (data.behavior.response > 0), gr);
% 
% if selectTrialNum
%     % only use those with a minimum trial count
%     minFun                = @(x, y) min([sum(x) sum(~x)]);
%     trialnum              = splitapply(minFun, (data.behavior.response > 0), gr);
%     minTrialCount         = 25;
%     biasedkernels(trialnum < minTrialCount, :) = NaN;
%     
%     scatter(1:length(trialnum), trialnum, 10, coh)
%     cnt = 1;
%     for c = unique(abs(coh))'
%         hold on
%         h{cnt} = histogram(trialnum(abs(coh) == c), 20, 'displaystyle', 'stairs');
%         cnt = cnt + 1;
%     end
%     legend([h{:}], {'0', '3', '9', '27', '81'});
% end
% 
% % average within subjects
% [gr, sj, pref]   = findgroups(sj, pref);
% meanFun          = @(x) nanmean(x, 1); % make sure to average over first dim
% biasedkernels    = splitapply(meanFun, biasedkernels, gr);
% 
% % then average over coherence levels within each subject
% subplot(4,4,[1 2]);
% hold on;
% colors = cbrewer('qual', 'Set1', 9);
% colors = colors([9 4], :);
% cohs = unique(pref);
% for c = 1:length(cohs),
%     plot(data.timeaxis, nanmean(biasedkernels(pref == cohs(c), :)), 'color', colors(c, :), 'linewidth', 0.5);
%     b{c} = boundedline(data.timeaxis(13:end), nanmean(biasedkernels(pref == cohs(c), 13:end)), ...
%         nanstd(biasedkernels(pref == cohs(c), 13:end)) ./ sqrt(length(unique(sj))), 'cmap', colors(c, :), 'alpha');
%     plot(data.timeaxis(13:end), nanmean(biasedkernels(pref == cohs(c), 13:end)), 'color', colors(c, :), 'linewidth', 1);
%     
% end
% legend([b{:}], {'non-preferred', 'preferred'}, 'location', 'eastoutside');
% legend boxoff;
% axis tight;
% 
% 
% %% do statistics on the timecourse
% % [h, p, stat] = ttest_clustercorr(biasedkernels(pref == 0, :), biasedkernels(pref == 1, :));
% [h, pval] = ttest(biasedkernels(pref == 0, :), biasedkernels(pref == 1, :));
% % [h, crit_p] = fdr_bh(pval, 0.05);
% 
% % remove significance during filter rise time
% h(1:12) = 0;
% 
% ylims = get(gca, 'Ylim');
% mask = double(h);
% mask(mask==0) = nan;
% mask = ((ylims(2)*0.3)+ylims(1))*mask; % plot a tiny bit above the lower ylim
% plot(data.timeaxis, mask, '.', 'MarkerSize', 10, 'color', 'k');
% 
% ylabel({'Excess motion'; 'energy fluctuations (%)'});
% xlabel('Time from stimulus onset (s)');
% axis tight; xlim([0 0.75]); set(gca, 'xtick', 0:0.25:0.75);
% box off; offsetAxes;
% 
% set(gca, 'xcolor', 'k', 'ycolor', 'k');
% tightfig;
% print(gcf, '-dpdf', '~/Data/serialHDDM/psychophysicalKernels_biased.pdf');

end