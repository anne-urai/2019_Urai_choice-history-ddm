function e1b_serialBias_SfN_ModelFreeCorrelation
  % from the csv table, make an overview of repetition behaviour

  % get a huge list with values for each participant
  % can then work with this dataframe

  clear; close all; clc;
  addpath(genpath('~/code/Tools'));
  warning off;

  % ============================================ %
  % TWO DIFFERENT DATASETS
  % ============================================ %

  plots = {'neutral', 'biased'};
  types = {'regress', 'stimcoding'};

  set(groot, 'defaultaxesfontsize', 6, 'defaultaxestitlefontsizemultiplier', 1, ...
  'defaultaxestitlefontweight', 'bold', ...
  'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters', ...
  'DefaultAxesBox', 'off', ...
  'DefaultAxesTickLength', [0.02 0.05], 'defaultaxestickdir', 'out', 'DefaultAxesTickDirMode', 'manual', ...
  'defaultfigurecolormap', [1 1 1], 'defaultTextInterpreter','tex');
  usr = getenv('USER');


  for s = 1:2,
    for p = 1:length(plots),

      switch p
      case 1

        switch usr
        case 'anne' % local
          datasets = {'RT_RDK', 'projects/0/neurodec/Data/MEG-PL', 'NatComm', 'Anke_2afc_neutral'};
        case 'aeurai' % lisa/cartesius
          datasets = {'NatComm', 'MEG', 'Anke_neutral', 'RT_RDK'};
        end
        datasetnames = {'2IFC (Urai et al. 2016)', '2IFC (MEG)', '2AFC (Braun et al.)', '2AFC (RT)'};

      case 2
        switch usr
        case 'aeurai' % lisa/cartesius
          datasets = {'Anke_alternating', 'Anke_neutral', 'Anke_repetitive'};
        end
        datasetnames = {'2AFC alternating', '2AFC neutral', '2AFC repetitive'};
      end

      % ============================================ %
      % ONE LARGE PLOT WITH PANEL FOR EACH DATASET
      % ============================================ %

      close all;
      for d = 1:length(datasets),

        if p == 2,
          colors = linspecer(3); % red blue green
          switch datasets{d}
          case 'Anke_alternating'
            plotColor = colors(1, :);
          case 'Anke_repetitive',
            plotColor = colors(3, :);
          case 'Anke_neutral'
            plotColor = colors(2, :);
          end
        else
          plotColor = [0.4 0.4 0.4];
        end

        results = readtable(sprintf('~/Data/%s/HDDM/summary/allindividualresults.csv', datasets{d}));

        % use the stimcoding difference
        switch types{s}
        case 'stimcoding'
          results.z_prevresp__regressdczprevresp = ...
          results.z_1__stimcodingdczprevresp - results.z_2__stimcodingdczprevresp;

          results.v_prevresp__regressdczprevresp = ...
          results.dc_1__stimcodingdczprevresp - results.dc_2__stimcodingdczprevresp;

        end

        subplot(4,4,d); hold on;
        plotScatter(results.v_prevresp__regressdczprevresp, results.repetition, 0.57, plotColor);
        title(datasetnames{d});
        xlabel('dc ~ prevresp');

        subplot(4,4,d+4); hold on;
        plotScatter(-1*results.z_prevresp__regressdczprevresp, results.repetition, 0.2, plotColor);
        xlabel('z ~ prevresp');

      end

      print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/figure1c_HDDM_modelfree_%s_%s.pdf', plots{p}, types{s}));
      print(gcf, '-depsc', sprintf('~/Data/serialHDDM/figure1c_HDDM_modelfree_%s_%s.eps', plots{p}, types{s}));
    end
  end
end

function plotScatter(x,y, legendWhere, plotColor);

  plot(x,y, '.');
  axis square;
  % axis tight;
  % % ylim([0.39 0.67]);
  % ticks = 0:0.1:1;
  % ticks(ticks < min(get(gca, 'ylim'))) = NaN;
  % ticks(ticks > max(get(gca, 'ylim'))) = NaN;
  % set(gca, 'ytick', ticks(~isnan(ticks)));
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
  offsetAxes;

  scatter(x,y, 10, ...
  'LineWidth', 0.01, ...
  'markeredgecolor', 'w', 'markerfacecolor', plotColor);
  ylabel('P(repeat)');

  txt = {sprintf('r = %.3f', rho) sprintf('p = %.3f', pval)};
  if pval < 0.001,
    txt = {sprintf('r = %.3f', rho) sprintf('p < 0.001')};
  end
  text(min(get(gca, 'xlim')) + legendWhere*(range(get(gca, 'xlim'))), ...
  min(get(gca, 'ylim')) + 0.15*(range(get(gca, 'ylim'))), ...
  txt, 'fontsize', 5);
end
