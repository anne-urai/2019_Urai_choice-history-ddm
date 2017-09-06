function e1b_serialBias_SfN_ModelFreeCorrelation
% from the csv table, make an overview of repetition behaviour

% get a huge list with values for each participant
% can then work with this dataframe

clear; close all; clc;
addpath(genpath('~/code/Tools'));

global mypath datasets datasetnames

% ============================================ %
% ONE LARGE PLOT WITH PANEL FOR EACH DATASET
% ============================================ %

close all;
for d = length(datasets):-1:1
    disp(datasets{d});

    colors = linspecer(5); % red blue green
    results = readtable(sprintf('%s/summary/%s/allindividualresults.csv', mypath, datasets{d}));
    results = results(results.session == 0, :);

    allresults = {};
    switch datasets{d}
        case {'Anke_2afc_sequential'}

            % DOES THIS DATASET HAVE MULTIPLE TRANSITION PROBABILITIES?
            % THEN PLOT THESE SEPARATELY

            % use the stimcoding difference only from alternating
            allresults{1}.z_prevresp        = results.z_1_neu__stimcodingdczprevresp - results.z_2_neu__stimcodingdczprevresp;
            allresults{1}.v_prevresp        = results.dc_1_neu__stimcodingdczprevresp - results.dc_2_neu__stimcodingdczprevresp;
            allresults{1}.criterionshift    = results.repetition_neutral;
            allresults{1}.subjnr            = results.subjnr;
            alltitles{1}                    = cat(2, datasetnames{d}{1}, ' - ', 'Neutral');

            allresults{2}.z_prevresp        = results.z_1_alt__stimcodingdczprevresp - results.z_2_alt__stimcodingdczprevresp;
            allresults{2}.v_prevresp        = results.dc_1_alt__stimcodingdczprevresp - results.dc_2_alt__stimcodingdczprevresp;
            allresults{2}.criterionshift    = results.repetition_alternating;
            allresults{2}.subjnr            = results.subjnr;
            alltitles{2}                    = cat(2, datasetnames{d}{1}, ' - ', 'Alternating');

            allresults{3}.z_prevresp        = results.z_1_rep__stimcodingdczprevresp - results.z_2_rep__stimcodingdczprevresp;
            allresults{3}.v_prevresp        = results.dc_1_neu__stimcodingdczprevresp - results.dc_2_rep__stimcodingdczprevresp;
            allresults{3}.criterionshift    = results.repetition_repetitive;
            allresults{3}.subjnr            = results.subjnr;
            alltitles{3}                    = cat(2, datasetnames{d}{1}, ' - ', 'Repetitive');

        case {'Bharath_fMRI', 'Anke_MEG', 'Anke_merged'}

            % DOES THIS DATASET HAVE MULTIPLE TRANSITION PROBABILITIES?
            % THEN PLOT THESE SEPARATELY

            % use the stimcoding difference only from alternating
            allresults{1}.z_prevresp        = results.z_1_0_50_0__stimcodingdczprevresp - results.z_2_0_50_0__stimcodingdczprevresp;
            allresults{1}.v_prevresp        = results.dc_1_0_50_0__stimcodingdczprevresp - results.dc_2_0_50_0__stimcodingdczprevresp;
            allresults{1}.criterionshift    = results.repetition_neutral;
            allresults{1}.subjnr            = results.subjnr;
            alltitles{1}                    = cat(2, datasetnames{d}{1}, ' - ', 'Neutral');

            allresults{2}.z_prevresp        = results.z_1_0_20_0__stimcodingdczprevresp - results.z_2_0_20_0__stimcodingdczprevresp;
            allresults{2}.v_prevresp        = results.dc_1_0_20_0__stimcodingdczprevresp - results.dc_2_0_20_0__stimcodingdczprevresp;
            allresults{2}.criterionshift    = results.repetition_alternating;
            allresults{2}.subjnr            = results.subjnr;
            alltitles{2}                    = cat(2, datasetnames{d}{1}, ' - ', 'Alternating');

            allresults{3}.z_prevresp        = results.z_1_0_80_0__stimcodingdczprevresp - results.z_2_0_80_0__stimcodingdczprevresp;
            allresults{3}.v_prevresp        = results.dc_1_0_80_0__stimcodingdczprevresp - results.dc_2_0_80_0__stimcodingdczprevresp;
            allresults{3}.criterionshift    = results.repetition_repetitive;
            allresults{3}.subjnr            = results.subjnr;
            alltitles{3}                    = cat(2, datasetnames{d}{1}, ' - ', 'Repetitive');

        otherwise

            try
                % use the stimcoding difference
                results.z_prevresp = ...
                    results.z_1__stimcodingdczprevresp - results.z_2__stimcodingdczprevresp;
                results.v_prevresp = ...
                    results.dc_1__stimcodingdczprevresp - results.dc_2__stimcodingdczprevresp;
            catch
                results.z_prevresp = ...
                    results.z_1_0__stimcodingdczprevresp - results.z_2_0__stimcodingdczprevresp;
                results.v_prevresp = ...
                    results.dc_1_0__stimcodingdczprevresp - results.dc_2_0__stimcodingdczprevresp;
            end

            results.criterionshift = results.repetition;
            allresults{1} = results;

            alltitles{1} = datasetnames{d}{1}; % use only the dataset title
    end


    for resultsSet = 1:length(allresults),

        results = allresults{resultsSet};

        % dont show errorbars for now
        if ~isfield(results, 'z_prevresp_ci'), results.z_prevresp_ci = nan(length(results.z_prevresp), 2); end
        if ~isfield(results, 'v_prevresp_ci'), results.v_prevresp_ci = nan(length(results.z_prevresp), 2); end
        if ~isfield(results, 'criterionshift_prct_1'), results.criterionshift_prct_1 = nan(length(results.z_prevresp), 1); end
        if ~isfield(results, 'criterionshift_prct_2'), results.criterionshift_prct_2 = nan(length(results.z_prevresp), 1); end

        disp(datasets{d}); disp(numel(unique(results.subjnr)));

        close all;
        subplot(4,4,1); hold on;
        [rho1, tt1] = plotScatter(results.z_prevresp, results.z_prevresp_ci, ...
            results.criterionshift, [results.criterionshift_prct_1 results.criterionshift_prct_2], 0.585, colors(5, :));
        xlabel('History bias in z', 'interpreter', 'tex');
        offsetAxes;
        if d == 2,
            set(gca, 'xtick', get(gca, 'xtick'), 'xticklabel', get(gca, 'xtick'));
        elseif d == 3 % avoid 10^-3 in the axis
            %   set(gca, 'xtick', get(gca, 'xtick'), 'xticklabel', {'', '0', '', '', '0.015'});
        end

        ylabel('\Deltac');
        ylabel('P(repeat)');

        sp2 = subplot(4,4,2); hold on;
        [rho2, tt2] = plotScatter(results.v_prevresp, results.v_prevresp_ci, ...
            results.criterionshift, [results.criterionshift_prct_1 results.criterionshift_prct_2], 0.05, colors(4, :));
        xlabel('History bias in v', 'interpreter', 'tex', 'fontweight', 'normal');
        set(gca, 'yticklabel', []);

        % compute the difference in correlation
        [rho3, pval3] = corr(results.v_prevresp, results.z_prevresp, ...
            'rows', 'complete', 'type', 'pearson');
        if pval3 < 0.05,
            fprintf('warning %s: rho = %.3f, pval = %.3f \n', datasets{d}, rho3, pval3);
        end
        [rhodiff, ~, pval] = rddiffci(rho1,rho2,rho3,numel(~isnan(results.criterionshift)), 0.05);
        offsetAxes; drawnow;

        % move together
        sp2.Position(1) = sp2.Position(1) - 0.08;
        ss = suplabel(alltitles{resultsSet}, 't');

        set(ss, 'fontweight', 'normal');
        ss.FontWeight = 'normal';
        ss.Position(2) = ss.Position(2) - 0.007;

        %% add line between the two correlation coefficients
        txt = {sprintf('\\Deltar_{%d} = %.3f, p = %.3f', length(find(~isnan(results.v_prevresp)))-3, rhodiff, pval)};
        if pval < 0.001,
            txt = {sprintf('\\Deltar_{%d} = %.3f, p < 0.001', length(find(~isnan(results.v_prevresp)))-3,  rhodiff)};
        end
        title(txt, 'fontweight', 'bold', 'fontsize', 6, 'horizontalalignment', 'left');
        tightfig;

        print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/figure1c_HDDM_modelfree_stimcoding_d%d_s%d.pdf', d, resultsSet));
    end
end

end

function [rho, tt] = plotScatter(x, xci, y, yci, legendWhere, plotColor);

plot(x,y, '.');
axis square;
axisNotSoTight;
[rho, pval] = corr(x,y, 'type', 'pearson', 'rows', 'complete');
l = lsline;

l.Color = 'k';
l.LineWidth = 0.5;
if pval < 0.05,
    l.LineStyle = '-';
else
    l.LineStyle = ':';
end

% show lines to indicate origin
xlims = [min(get(gca, 'xlim')) max(get(gca, 'xlim'))];
ylims = [min(get(gca, 'ylim')) max(get(gca, 'ylim'))];
plot([0 0], ylims, 'color', [0.5 0.5 0.5], 'linewidth', 0.2);
plot(xlims, [0.5 0.5], 'color', [0.5 0.5 0.5], 'linewidth', 0.2); % if p(repeat), 0.5

% % scatter with error bars
% h = ploterr(x, y, xci-x, yci-y, '.', 'abshhxy', 0);
%
% % color markers appropiately
% set([h(1)], 'markersize', 8, ...
%     'LineWidth', 0.001, ...
%     'markeredgecolor', 'w', 'markerfacecolor', plotColor);
% set([h(2:end)], 'color', plotColor, 'linewidth', 0.2);

s  = scatter(x, y,  15, ...
    'LineWidth', 0.001, ...
    'markeredgecolor', 'w', 'markerfacecolor', plotColor);

txt = {sprintf('r_{%d} = %.3f', length(find(~isnan(y)))-2, rho) sprintf('p = %.3f', pval)};
if pval < 0.001,
    txt = {sprintf('r_{%d} = %.3f', length(find(~isnan(y)))-2,rho) sprintf('p < 0.001')};
end
tt = text(min(get(gca, 'xlim')) + legendWhere*(range(get(gca, 'xlim'))), ...
    min(get(gca, 'ylim')) + 0.8*(range(get(gca, 'ylim'))), ...
    txt, 'fontsize', 5);
set(gca, 'color', 'none');

% also add the group mean
p = ploterr(nanmean(x), nanmean(y), nanstd(x) ./ sqrt(length(x)), ...
    nanstd(y) ./ sqrt(length(y)), '.k', 'abshhxy', 0);
set(p(1), 'markersize', 0.1); % tiny marker

end
