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

for maxcohlevel = [0,3,9,27,81],
    
    close all;
    
    % either all the data or only neutral
    load(sprintf('%s/%s', path, 'motionEnergyData_AnkeMEG.mat'));
    
    % use the normalized motion energy, so that the units are % coherence 'up'
    data.motionenergy = data.motionenergy_normalized;
    
    % only select trials without any objective evidence level - remove high
    % coherence trials!
    data.motionenergy = data.motionenergy([data.behavior.coherence] <= maxcohlevel, :);
    data.behavior     = data.behavior([data.behavior.coherence] <= maxcohlevel, :);
    
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
    kernels               = splitapply(kernelFun, data.motionenergy, (data.behavior.response > 0), gr);
    
    % average within each subject!
    if numel(unique(coh)) > 1,
        kernels = splitapply(@nanmean, kernels, findgroups(sj));
    end
    
    % =============================== %
    % PSYCHOPHYSICAL KERNELS - without bias
    % =============================== %
    
    % then average over coherence levels within each subject
    subplot(441);
    hold on;
    colors(1, :) = [0.3 0.3 0.3]; c = 1;
    plot(data.timeaxis, nanmean(kernels), 'color', colors(c, :), 'linewidth', 0.5);
    b{c} = boundedline(data.timeaxis(13:end), nanmean(kernels(:, 13:end)), ...
        nanstd(kernels(:, 13:end)) ./ sqrt(length(unique(sj))), 'cmap', colors(c, :), 'alpha');
    plot(data.timeaxis(13:end), nanmean(kernels(:, 13:end)), 'color', colors(c, :), 'linewidth', 1);
    axis tight;
    
    % do statistics on the timecourse
    [h, p, stat] = ttest_clustercorr(kernels(:, 13:end));
    %[h, pval] = ttest(biasedkernels);
    %[h, crit_p] = fdr_bh(pval, 0.05);
    
    ylims = get(gca, 'Ylim');
    mask = [zeros(1,12) double(h)];
    mask(mask==0) = nan;
    mask = ((ylims(2)*0.1)+ylims(1))*mask; % plot a tiny bit above the lower ylim
    plot(data.timeaxis, mask, '.', 'MarkerSize', 10, 'color', 'k');
    
    ylabel({'Excess motion'; 'energy fluctuations (%)'});
    xlabel('Time from stimulus onset (s)');
    axis tight; xlim([0 0.75]); set(gca, 'xtick', 0:0.25:0.75);
    box off; ylim([-0.5 2]); offsetAxes;
    set(gca, 'xcolor', 'k', 'ycolor', 'k');
    %tightfig;
    %print(gcf, '-dpdf', '~/Data/serialHDDM/psychophysicalKernels.pdf');
    
    
    % %% COMPARE THE KERNELS WITH THE O-U EFFECTIVE LEAK PARAMETER
    %
    % global mypath;
    % results     = readtable(sprintf('%s/summary/%s/allindividualresults_kostis.csv', mypath, 'Anke_MEG_Neutral'));
    % lambda      = results.ouK_vanilla_dv;
    % kernelDiff  = nanmean(kernels(:, 13:29), 2) - nanmean(kernels(:, 30:45), 2);
    % [rho, pval] = corr(kernelDiff, lambda, 'type', 'spearman');
    % fprintf('\n\nCorrelation between kernels and O-U lambda, Spearmans rho: %.3f, p = %.3f \n', rho, pval);
    
    %% ============================================================== %
    % PSYCHOPHYSICAL KERNELS - repetition vs. alternation, for 2 subgroups
    % ============================================================== %
    
    % bin
    coh                   = [data.behavior.stimulus .* data.behavior.coherence];
    [gr, sj, rep, coh]    = findgroups(data.behavior.subj_idx, data.behavior.repeat, coh);
    kernels               = splitapply(kernelFun, data.motionenergy, (data.behavior.response > 0), gr);
    
    % average within each subject and bias
    [gr, sj, rep] = findgroups(sj, rep);
    if numel(unique(coh)) > 1,
        kernels = splitapply(@nanmean, kernels, gr);
    end
    
    kernels_rep = kernels(rep == 1, :);
    kernels_alt = kernels(rep == 0, :);
    
    sjrep       = splitapply(@nanmean, data.behavior.repeat, findgroups(data.behavior.subj_idx));
    group.alternators = (sjrep < 0.5);
    group.repeaters   = (sjrep > 0.5);
    groups = {'repeaters', 'alternators'};
    
    for g = 1:length(groups),
        % then average over coherence levels within each subject
        subplot(4,4,g+1);
        hold on;
        
        colors = cbrewer('qual', 'Set1', 5);
        colors = colors([4 5], :);
        
        % repetition kernels
        plot(data.timeaxis, nanmean(kernels_rep(group.(groups{g}), :)), 'color', colors(1, :), 'linewidth', 0.5);
        b{1} = boundedline(data.timeaxis(13:end), nanmean(kernels_rep(group.(groups{g}), 13:end)), ...
            nanstd(kernels_rep(group.(groups{g}), 13:end)) ./ sqrt(sum((group.(groups{g})))), 'cmap', colors(1, :), 'alpha');
        plot(data.timeaxis(13:end), nanmean(kernels_rep(group.(groups{g}), 13:end)), 'color', colors(1, :), 'linewidth', 1);
        
        % alternation kernels
        b{2} = plot(data.timeaxis, nanmean(kernels_alt(group.(groups{g}), :)), 'color', colors(2, :), 'linewidth', 0.5);
        boundedline(data.timeaxis(13:end), nanmean(kernels_alt(group.(groups{g}), 13:end)), ...
            nanstd(kernels_alt(group.(groups{g}), 13:end)) ./ sqrt(sum((group.(groups{g})))), 'cmap', colors(2, :), 'alpha');
        plot(data.timeaxis(13:end), nanmean(kernels_alt(group.(groups{g}), 13:end)), 'color', colors(2, :), 'linewidth', 1);
        
        % do statistics on the timecourse
        [h, p, stat] = ttest_clustercorr(kernels_alt(group.(groups{g}), 13:end), kernels_rep(group.(groups{g}), 13:end));
        
        axis tight;
        ylims = get(gca, 'Ylim');
        mask = [zeros(1,12) double(h)];
        mask(mask==0) = nan;
        mask = ((ylims(2)*0.1)+ylims(1))*mask; % plot a tiny bit above the lower ylim
        plot(data.timeaxis, mask, '.', 'MarkerSize', 10, 'color', 'k');
        
        % if g == 1, ylabel({'Excess motion'; 'energy fluctuations (%)'}); end
        xlabel('Time from stimulus onset (s)');
        axis tight; xlim([0 0.75]); set(gca, 'xtick', 0:0.25:0.75);
        box off; % ylim([-0.5 2]); 
        offsetAxes;
        set(gca, 'xcolor', 'k', 'ycolor', 'k');
        title(capitalize(groups{g}));
    end
    
    l = legend([b{:}], {'repeat', 'alternate'});
    l.Position(2) = l.Position(2) - 0.2;
    
    %     sp = subplot(4,8,5);
    %     axis off;
    %     l.Position = get(sp, 'Position');
    l.Box = 'off';
    
    % tightfig;
    % print(gcf, '-dpdf', '~/Data/serialHDDM/psychophysicalKernels_repeat_alternate.pdf');
    
    
    %% ============================================================== %
    % PSYCHOPHYSICAL KERNELS - collapse repeat and alternate into bias
    % ============================================================== %
    
    kernels_pref = [kernels_rep(group.repeaters, :); kernels_alt(group.alternators, :)];
    kernels_unpref = [kernels_alt(group.repeaters, :); kernels_rep(group.alternators, :)];
    
    % then average over coherence levels within each subject
    %close all;
    subplot(4,4,4);
    hold on;
    
    colors = cbrewer('qual', 'Set1', 9);
    colors = colors([8 9], :);
    
    % repetition kernels
    plot(data.timeaxis, nanmean(kernels_pref), 'color', colors(1, :), 'linewidth', 0.5);
    b{1} = boundedline(data.timeaxis(13:end), nanmean(kernels_pref(:, 13:end)), ...
        nanstd(kernels_pref(:, 13:end)) ./ sqrt(numel(unique(data.behavior.subj_idx))), 'cmap', colors(1, :), 'alpha');
    plot(data.timeaxis(13:end), nanmean(kernels_pref(:, 13:end)), 'color', colors(1, :), 'linewidth', 1);
    
    % alternation kernels
    plot(data.timeaxis, nanmean(kernels_unpref), 'color', colors(2, :), 'linewidth', 0.5);
    b{2} = boundedline(data.timeaxis(13:end), nanmean(kernels_unpref(:, 13:end)), ...
        nanstd(kernels_unpref(:, 13:end)) ./ sqrt(numel(unique(data.behavior.subj_idx))), 'cmap', colors(2, :), 'alpha');
    plot(data.timeaxis(13:end), nanmean(kernels_unpref(:, 13:end)), 'color', colors(2, :), 'linewidth', 1);
    
    % do statistics on the timecourse
    [h, p, stat] = ttest_clustercorr(kernels_unpref(:, 13:end), kernels_pref(:, 13:end));
    
    axis tight;
    ylims = get(gca, 'Ylim');
    mask = [zeros(1,12) double(h)];
    mask(mask==0) = nan;
    mask = ((ylims(2)*0.1)+ylims(1))*mask; % plot a tiny bit above the lower ylim
    plot(data.timeaxis, mask, '.', 'MarkerSize', 10, 'color', 'k');
    
    % ylabel({'Excess motion'; 'energy fluctuations (%)'});
    xlabel('Time from stimulus onset (s)');
    axis tight; xlim([0 0.75]); set(gca, 'xtick', 0:0.25:0.75);
    box off; ylim([-0.5 2]); offsetAxes;
    set(gca, 'xcolor', 'k', 'ycolor', 'k');
    
    l = legend([b{:}], {'preferred', 'unpreferred'});
    l.Position(2) = l.Position(2) - 0.2;
    %     sp = subplot(4,8,5);
    %     axis off;
    %     l.Position = get(sp, 'Position');
    l.Box = 'off';
    
    %tightfig;
    % print(gcf, '-dpdf', '~/Data/serialHDDM/psychophysicalKernels_pref_unpref.pdf');
    
    suplabel(sprintf('Coherence levels up to %d%% included', maxcohlevel), 't');
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/psychophysicalKernels_maxCoh%d.pdf', maxcohlevel));
end

end