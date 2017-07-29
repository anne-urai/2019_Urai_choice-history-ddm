function e1b_serialBias_SfN_ModelFreeCorrelation
% from the csv table, make an overview of repetition behaviour

% get a huge list with values for each participant
% can then work with this dataframe

clear; close all; clc;
addpath(genpath('~/code/Tools'));

global mypath datasets datasetnames

% ============================================ %
% TWO DIFFERENT DATASETS
% ============================================ %

% ============================================ %
% ONE LARGE PLOT WITH PANEL FOR EACH DATASET
% ============================================ %

close all;
for d = 1:length(datasets),

    colors = linspecer(5); % red blue green

    results = readtable(sprintf('%s/%s/summary/allindividualresults.csv', mypath, datasets{d}));
    results = results(results.session == 0, :);
    disp(datasets{d}); disp(numel(unique(results.subjnr)));

    % use the stimcoding difference

    results.z_prevresp__regressdczprevresp = ...
        results.z_1__stimcodingdczprevresp - results.z_2__stimcodingdczprevresp;
    results.v_prevresp__regressdczprevresp = ...
        results.dc_1__stimcodingdczprevresp - results.dc_2__stimcodingdczprevresp;

    close all;
    subplot(4,4,1); hold on;
    rho1 = plotScatter(results.v_prevresp__regressdczprevresp, results.repetition, 0.57, colors(4, :));
    title(datasetnames{d});
    ll = xlabel('dc ~ previous response');
    offsetAxes;

     switch d
        case {1, 5}
            ylabel('P(repeat)');
        case 2
            set(gca, 'ytick', 0.46:0.04:0.58);
    end

    sp2 = subplot(4,4,5); hold on;
    rho2 = plotScatter(results.z_prevresp__regressdczprevresp, results.repetition, 0.57, colors(5, :));
    xlabel('z ~ previous response');

    switch d
        case {1, 5}
            ylabel('P(repeat)');
        case 2
            set(gca, 'ytick', 0.46:0.04:0.58);
    end

    % compute the difference in correlation
    rho3 = corr(results.v_prevresp__regressdczprevresp, results.z_prevresp__regressdczprevresp, ...
        'rows', 'complete', 'type', 'pearson');
    [rhodiff, cihilow, pval] = rddiffci(rho1,rho2,rho3,numel(~isnan(results.repetition)), 0.05);
    % disp(rho3);

    txt = sprintf('\\Deltar = %.3f, p = %.3f', rhodiff, pval);
    if pval < 0.001,
        txt = sprintf('\\Deltar = %.3f, p < 0.001', rhodiff);
    end
    title({' '; txt}, 'fontweight', 'normal', 'fontsize', 5);
    offsetAxes; drawnow;
  %  sp2.Position(2) = sp2.Position(2) - 0.01;

    tightfig;
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/figure1c_HDDM_modelfree_stimcoding_d%d.pdf',d));
end
close all;

end

function rho = plotScatter(x,y, legendWhere, plotColor);

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
plot(xlims, [0.5 0.5], 'color', [0.5 0.5 0.5], 'linewidth', 0.2);

scatter(x,y, 15, ...
    'LineWidth', 0.001, ...
    'markeredgecolor', 'w', 'markerfacecolor', plotColor);

txt = {sprintf('r = %.3f', rho) sprintf('p = %.3f', pval)};
if pval < 0.001,
    txt = {sprintf('r = %.3f', rho) sprintf('p < 0.001')};
end
text(min(get(gca, 'xlim')) + legendWhere*(range(get(gca, 'xlim'))), ...
    min(get(gca, 'ylim')) + 0.15*(range(get(gca, 'ylim'))), ...
    txt, 'fontsize', 5);
    set(gca, 'color', 'none');

end
