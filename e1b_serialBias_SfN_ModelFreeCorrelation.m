% from the csv table, make an overview of repetition behaviour

% get a huge list with values for each participant
% can then work with this dataframe

clear; close all; clc;
addpath(genpath('~/code/Tools'));
warning off;

% ============================================ %
% TWO DIFFERENT DATASETS
% ============================================ %

usr = getenv('USER');
switch usr
case 'anne' % local
  datasets = {'RT_RDK', 'projects/0/neurodec/Data/MEG-PL', 'NatComm', 'Anke_2afc_neutral'};
case 'aeurai' % lisa/cartesius
  datasets = {'NatComm', 'MEG', 'Anke_neutral', 'RT_RDK'};
end
datasetnames = {'2IFC (Urai et al. 2016)', '2IFC (MEG)', '2AFC (Braun et al.)', '2AFC (RT)'};

set(groot, 'defaultaxesfontsize', 6, 'defaultaxestitlefontsizemultiplier', 1, ...
    'defaultaxestitlefontweight', 'bold', ...
    'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'opengl', ...
    'DefaultAxesBox', 'off', ...
    'DefaultAxesTickLength', [0.02 0.05], 'defaultaxestickdir', 'out', 'DefaultAxesTickDirMode', 'manual', ...
    'defaultfigurecolormap', [1 1 1]);

if 0,
for d = 1:length(datasets),
    results = readtable(sprintf('~/Data/%s/HDDM/summary/allindividualresults.csv', datasets{d}));

    % ============================================ %
    % compute repetition parameters from separate HDDM models
    % ============================================ %

    results.dc_prevresp__stimcodingdczprevresp = ...
        results.dc_1__stimcodingdczprevresp - results.dc_2__stimcodingdczprevresp;

    results.z_prevresp__stimcodingdczprevresp = ...
        results.z_1__stimcodingdczprevresp - results.z_2__stimcodingdczprevresp;

    % ============================================ %
    % RENAME PARAMETERS
    % ============================================ %

    results.Properties.VariableNames{'dc_prevresp__stimcodingdczprevresp'}      = 'dc_seq_stimcoding';
    results.Properties.VariableNames{'v_prevresp__regressdczprevresp'}          = 'dc_seq_regress';
    results.Properties.VariableNames{'z_prevresp__stimcodingdczprevresp'}       = 'z_seq_stimcoding';
    results.Properties.VariableNames{'z_prevresp__regressdczprevresp'}          = 'z_seq_regress';

    % ============================================ %
    % SEPARATE OR JOINT FIT
    % ============================================ %

    close;
    corrplot(results, {'dc_seq_regress', 'z_seq_regress', 'dc_seq_stimcoding', 'z_seq_stimcoding'}, ...
        {'repetition'});
    suplabel(sprintf('%s', datasetnames{d}), 't');
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/%s_corrplot.pdf', datasetnames{d}));

end
end

% ============================================ %
% ONE LARGE PLOT WITH PANEL FOR EACH DATASET
% ============================================ %

close all;
for d = 1:length(datasets),
    results = readtable(sprintf('~/Data/%s/HDDM/summary/allindividualresults.csv', datasets{d}));

    % ============================================ %
    % DC
    % ============================================ %

    subplot(4,4,d); hold on;
    plot(results.v_prevresp__regressdczprevresp, results.repetition, '.');
    axis square; axisNotSoTight;
    ylim([0.39 0.67]);
    set(gca, 'ytick', [0.4 0.5 0.6]);

    [rho, pval] = corr(results.v_prevresp__regressdczprevresp, results.repetition, 'type', 'pearson', 'rows', 'complete');
    l = lsline;
    l.Color = 'k';
    l.LineWidth = 0.5;
    if pval < 0.05,
      l.LineStyle = '-';
    else
      l.LineStyle = ':';
    end

    scatter(results.v_prevresp__regressdczprevresp, results.repetition, 15, ...
    'LineWidth', 0.01, ...
    'markeredgecolor', 'w', 'markerfacecolor', [0.4 0.4 0.4]);
    title(datasetnames{d});
    ylabel('P(repeat)');
    xlabel('dc ~ resp_{-1}');
    offsetAxes;

    txt = {sprintf('r = %.3f', rho) sprintf('p = %.3f', pval)};
    if pval < 0.001,
      txt = {sprintf('r = %.3f', rho) sprintf('p < 0.001')};
    end
    text(min(get(gca, 'xlim')) + 0.57*(range(get(gca, 'xlim'))), ...
    min(get(gca, 'ylim')) + 0.15*(range(get(gca, 'ylim'))), ...
    txt, 'fontsize', 5);

    % add boxplot to show distribution

    % ============================================ %
    % Z
    % ============================================ %

    subplot(4,4,d+4); hold on;
    plot(results.z_prevresp__regressdczprevresp, results.repetition, '.');
    axis square; axisNotSoTight;
    ylim([0.39 0.67]);
    set(gca, 'ytick', [0.4 0.5 0.6]);
    [rho, pval] = corr(results.z_prevresp__regressdczprevresp, results.repetition, 'type', 'pearson', 'rows', 'complete');
    l = lsline;
    l.Color = 'k';
    l.LineWidth = 0.5;
    if pval < 0.05,
      l.LineStyle = '-';
    else
      l.LineStyle = ':';
    end

    scatter(results.z_prevresp__regressdczprevresp, results.repetition, 15, ...
    'LineWidth', 0.01, ...
    'markeredgecolor', 'w', 'markerfacecolor', [0.4 0.4 0.4]);
    % title(datasetnames{d});
    ylabel('P(repeat)');
    xlabel('z ~ resp_{-1}');
    offsetAxes;

    txt = {sprintf('r = %.3f', rho) sprintf('p = %.3f', pval)};
    if pval < 0.001,
      txt = {sprintf('r = %.3f', rho) sprintf('p < 0.001')};
    end
    text(min(get(gca, 'xlim')) + 0.57*(range(get(gca, 'xlim'))), ...
    min(get(gca, 'ylim')) + 0.15*(range(get(gca, 'ylim'))), ...
    txt, 'fontsize', 5);
    disp(get(gca, 'ylim'));
  end

  print(gcf, '-depsc', sprintf('~/Data/serialHDDM/HDDM_modelfree_correlation.eps'));
  %export_fig(sprintf('~/Data/serialHDDM/HDDM_modelfree_correlation.pdf'));
