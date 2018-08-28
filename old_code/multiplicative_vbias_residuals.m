function multiplicative_vbias_residuals

% Code to fit the history-dependent drift diffusion models described in
% Urai AE, Gee JW de, Donner TH (2018) Choice history biases subsequent evidence accumulation. bioRxiv:251595
%
% MIT License
% Copyright (c) Anne Urai, 2018
% anne.urai@gmail.com

close all; clc;
addpath(genpath('~/code/Tools'));

global mypath datasets datasetnames colors

% ============================================ %
% ONE LARGE PLOT WITH PANEL FOR EACH DATASET
% ============================================ %

for d = [4 5] % only NatComm and Anke_MEG_neutral, with varying coherence level
    disp(datasets{d});
    
    % load simulated data - make sure this has all the info we need
    alldata    = readtable(sprintf('%s/summary/%s/%s_ppc_data.csv', mypath, datasets{d}, 'stimcoding_dc_prevresp'));
    alldata    = sortrows(alldata, {'subj_idx'});
    
    % recode into repeat and alternate for the model
    alldata.biased = ((alldata.response > 0) == (alldata.prevresp > 0));
    alldata.biased_model = ((alldata.response_sampled > 0) == (alldata.prevresp > 0));
    
    % for each observers, compute their bias
    [gr, sjs] = findgroups(alldata.subj_idx);
    sjrep = splitapply(@nanmean, alldata.biased, gr);
    sjrep = sjs(sjrep < 0.5);
    
    % recode into biased and unbiased choices
    altIdx                      = ismember(alldata.subj_idx, sjrep);
    alldata.biased(altIdx)      = double(~(alldata.biased(altIdx))); % flip
    alldata.biased_model(altIdx) = double(~(alldata.biased_model(altIdx))); % flip
    
    % FOR BOTH REAL AND SIMULATED CHOICES, SEE BIAS PROPORTION AS A
    % FUNCTION OF COHERENCE
    [gr, sj, coh]   = findgroups(alldata.subj_idx, alldata.coherence);
    tab = array2table([sj, coh], 'variablenames', {'subj_idx', 'coherence'});
    tab_data = tab;     tab_model = tab;
    tab_data.bias   = splitapply(@nanmean, alldata.biased, gr);
    tab_model.bias  = splitapply(@nanmean, alldata.biased_model, gr);
    mat_data        = unstack(tab_data, 'bias', 'coherence');
    mat_model       = unstack(tab_model, 'bias', 'coherence');

    % PLOT
    subplot(3,3,d); hold on;
    cohlevels = unique(abs(alldata.coherence)) * 100;
    boundedline(cohlevels, ...
        nanmean(mat_data{:, 2:end}), nanstd(mat_data{:, 2:end}) ./ sqrt(size(mat_data, 1)), 'alpha');
    boundedline(cohlevels, ...
        nanmean(mat_model{:, 2:end}), nanstd(mat_model{:, 2:end}) ./ sqrt(size(mat_data, 1)), 'cmap', [0 0 0], 'alpha');
    set(gca, 'xtick', cohlevels, 'xticklabelrotation', -30);
    xlabel('Coherence (%)');
    ylabel('P(bias)');
    box off; offsetAxes;
    title(datasetnames{d});
    
    subplot(3,3,d+3); hold on;
    diff = mat_data{:, 2:end} - mat_model{:, 2:end};
    boundedline(cohlevels, ...
        nanmean(diff), nanstd(diff) ./ sqrt(size(mat_data, 1)), 'alpha');
    set(gca, 'xtick', cohlevels, 'xticklabelrotation', -30);
    xlabel('Coherence (%)');
    ylabel('Residuals');
    box off; offsetAxes;
    
    % SHOW STATS
    [h, pval] = ttest(mat_data{:, 2:end}, mat_model{:, 2:end});
    [h, crit_p] = fdr_bh(pval, 0.05);
    
    ylims = get(gca, 'Ylim');
    mask = double(h);
    mask(mask==0) = nan;
    mask = ((ylims(2)*0.3)+ylims(1))*mask; % plot a tiny bit above the lower ylim
    plot(cohlevels, mask, '.', 'MarkerSize', 10, 'color', 'k');
    
end


tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/vbias_per_coherence_residuals.pdf'));
