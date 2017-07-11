function e3_serialBias_SfN_RTmodulation

  addpath(genpath('~/code/Tools'));
  warning off; close all; clear;
  global datasets datasetnames

  set(groot, 'defaultaxesfontsize', 6, 'defaultaxestitlefontsizemultiplier', 1, ...
  'defaultaxestitlefontweight', 'bold', ...
  'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters', ...
  'DefaultAxesBox', 'off', ...
  'DefaultAxesTickLength', [0.02 0.05], 'defaultaxestickdir', 'out', 'DefaultAxesTickDirMode', 'manual', ...
  'defaultfigurecolormap', [1 1 1], 'defaultTextInterpreter','tex');
  usr = getenv('USER');

  % ========================================== %
  % MODULATION OF SERIAL CHOICE BIAS
  % ========================================== %

  plots = {'neutral', 'biased'};

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

    cnt = 1; close all;
    for d = 1:length(datasets),

      % get traces for the model with pupil and rt modulation
      traces = readtable(sprintf('~/Data/%s/HDDM/regress_dc_prevresp_prevstim_prevrt_prevpupil/all_traces.csv', datasets{d}));

      colors = cbrewer('qual', 'Paired', 8);
      % plot the pupil and RT traces
      subplot(4,4,d); hold on;
      h1 = histogram_smooth(traces.v_prevresp_prevrt, colors(1, :), colors(2, :));
      h2 = histogram_smooth(traces.v_prevresp_prevpupil, colors(3, :), colors(4, :));

      % show if these are significant (1-sided?)
      pvalRT = mean(traces.v_prevresp_prevrt > 0);
      pvalPupil = mean(traces.v_prevresp_prevpupil > 0);

      axis tight; axis square;
      l = legend([h1 h2], {sprintf('RT, p = %.3f', pvalRT); sprintf('Pupil, p = %.3f', pvalPupil)}, ...
      'location', 'southeast');
      l.Position(2) = l.Position(2) - 0.15;
      legend boxoff;
      title(datasetnames{d}); xlabel('dc ~ prevresp modulation');
      vline(0);
      offsetAxes_y;
    end

    print(gcf, '-depsc', sprintf('~/Data/serialHDDM/figure2_modulationTraces_%s.eps', plots{p}));
  end
end

function h = histogram_smooth(x, color1, color2)
  [f,xi] = ksdensity(x);
  area(xi, f, 'edgecolor', 'none', 'facecolor', color1);
  % area
  h = plot(xi, f, 'color', color2, 'linewidth', 1);
end

function offsetAxes_y()

  if ~exist('ax', 'var'), ax = gca;
  end
  if ~exist('offset', 'var'), offset = 4;
  end

  % ax.YLim(1) = ax.YLim(1)-(ax.YTick(2)-ax.YTick(1))/offset;
  ax.XLim(1) = ax.XLim(1)-(ax.XTick(2)-ax.XTick(1))/offset;

  % this will keep the changes constant even when resizing axes
  addlistener(ax, 'MarkedClean', @(obj,event)resetVertex(ax));
end

function resetVertex ( ax )
  % repeat for Y (set 2nd row)
  ax.YRuler.Axle.VertexData(2,1) = min(get(ax, 'Ytick'));
  ax.YRuler.Axle.VertexData(2,2) = max(get(ax, 'Ytick'));
  % X, Y and Z row of the start and end of the individual axle.
  ax.XRuler.Axle.VertexData(1,1) = min(get(ax, 'Xtick'));
  ax.XRuler.Axle.VertexData(1,2) = max(get(ax, 'Xtick'));
end
