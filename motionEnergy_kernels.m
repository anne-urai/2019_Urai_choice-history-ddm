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

path = '~/Data/psychophysicalKernels';
load(sprintf('%s/%s', path, 'motionEnergyData_AnkeMEG.mat'));

% only select neutral
% data.motionenergy = data.motionenergy([data.behavior.transitionprob] == 0.5, :);
% data.behavior     = data.behavior([data.behavior.transitionprob] == 0.5, :);

% recode into choices that are biased
data.behavior.repeat = (data.behavior.response == data.behavior.prevresp);

% for each observers, compute their bias
[gr, sjs]   = findgroups(data.behavior.subj_idx);
sjrep       = splitapply(@nanmean, data.behavior.repeat, gr);
sjrep       = sjs(sjrep < 0.5);

% recode into biased and unbiased choices
data.behavior.biased = data.behavior.repeat;
altIdx      = ismember(data.behavior.subj_idx, sjrep);
data.behavior.biased(altIdx) = double(~(data.behavior.biased(altIdx))); % flip

% =============================== %
% PSYCHOPHYSICAL KERNELS
% =============================== %

meanSubtract = 1;
% only select trials without any objective evidence
%if meanSubtract,
data.motionenergy = data.motionenergy_normalized - data.behavior.stimulus .* data.behavior.coherence;
%else
% only select trials without any objective evidence
data.motionenergy = data.motionenergy([data.behavior.coherence] < 10, :);
data.behavior     = data.behavior([data.behavior.coherence] < 10, :);
%end

% what is the time-course of evidence that leads subjects to make their
% preferred vs non-preferred choice?
kernelFun       = @(x, y) nanmean(x(y, :)) - nanmean(x(~y, :));

% average those within each subject?
close all;
plot(data.timeaxis, kernelFun(data.motionenergy, (data.behavior.response > 0)));

% bin
coh                   = [data.behavior.stimulus .* data.behavior.coherence];
[gr, sj, coh, pref]   = findgroups(data.behavior.subj_idx, coh, data.behavior.biased);
biasedkernels         = splitapply(kernelFun, data.motionenergy, (data.behavior.response > 0), gr);

[gr, sj, pref]   = findgroups(sj, pref);
biasedkernels    = splitapply(@nanmean, biasedkernels, gr);

% then average over coherence levels within each subject
subplot(4,4,[1 2]); 
hold on;
colors = cbrewer('qual', 'Set1', 9);
colors = colors([9 4], :);
cohs = unique(pref);
for c = 1:length(cohs),
    plot(data.timeaxis, nanmean(biasedkernels(pref == cohs(c), :)), 'color', colors(c, :), 'linewidth', 0.5);
    b{c} = boundedline(data.timeaxis(13:end), nanmean(biasedkernels(pref == cohs(c), 13:end)), ...
        nanstd(biasedkernels(pref == cohs(c), 13:end)) ./ sqrt(length(unique(sj))), 'cmap', colors(c, :), 'alpha');
    plot(data.timeaxis(13:end), nanmean(biasedkernels(pref == cohs(c), 13:end)), 'color', colors(c, :), 'linewidth', 1);

end
legend([b{:}], {'non-preferred', 'preferred'}, 'location', 'eastoutside');
legend boxoff;

%% do statistics on the timecourse
% [h, p, stat] = ttest_clustercorr(biasedkernels(pref == 0, :), biasedkernels(pref == 1, :));
[h, pval] = ttest(biasedkernels(pref == 0, :), biasedkernels(pref == 1, :));
% [h, crit_p] = fdr_bh(pval, 0.05);

ylims = get(gca, 'Ylim');
mask = double(h);
mask(mask==0) = nan;
mask = ((ylims(2)*0.3)+ylims(1))*mask; % plot a tiny bit above the lower ylim
plot(data.timeaxis, mask, '.', 'MarkerSize', 10, 'color', 'k');

% boundedline(data.timeaxis, nanmean(biasedkernels), ...
%     nanstd(biasedkernels) ./ sqrt(length(unique(sj))), 'cmap', [0 0 0], 'alpha');
% plot(data.timeaxis(13:end), nanmean(biasedkernels(:, 13:end)), 'k', 'linewidth', 2);
ylabel({'Excess motion'; 'energy fluctuations (%)'});
xlabel('Time from stimulus onset (s)');
axis tight; xlim([0 0.75]); set(gca, 'xtick', 0:0.25:0.75);
box off; offsetAxes;

set(gca, 'xcolor', 'k', 'ycolor', 'k');
tightfig;
print(gcf, '-dpdf', '~/Data/serialHDDM/psychophysicalKernels.pdf');

end