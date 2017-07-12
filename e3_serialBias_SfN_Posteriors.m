function e3_serialBias_SfN_Posteriors

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

  inverse_link_function = @(x) 1/(1+exp(-x));
  %inverse_link_function = @(x) -x; %1/(1-exp(-x));

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
    parameters = {'dc', 'z'}; cnt = 1;
    for pa = 1:length(parameters),
      for d = 1:length(datasets),

        % get traces for the model with pupil and rt modulation
        traces = readtable(sprintf('~/Data/%s/HDDM/regress_dc_z_prevresp_prevstim/all_traces.csv', datasets{d}));

        switch parameters{pa}
          case 'dc'
          dat1 = traces.v_prevresp;
          dat2 = traces.v_prevstim;
          case 'z'
          dat1 = -traces.z_prevresp;
          dat2 = -traces.z_prevstim;
        end

        % plot the pupil and RT traces
        purples = cbrewer('seq', 'Purples', 3);
        greys = cbrewer('seq', 'Greys', 3);

        subplot(4,4,cnt); hold on; cnt = cnt + 1;
        h1 = histogram_smooth(dat1, greys(1, :), greys(end, :));
        h2 = histogram_smooth(dat2, purples(1, :), purples(end, :));

        % show if these are significant - two sided
        % https://github.com/jwdegee/2017_eLife/blob/master/hddm_regression.py, line 273
        pval1    = min([mean(dat1 > 0) mean(dat1 < 0)]);
        pval2    = min([mean(dat2 > 0) mean(dat2 < 0)]);

        axis tight; axis square;
        l = legend([h1 h2], {sprintf('prevresp, p = %.3f', pval1); sprintf('prevstim, p = %.3f', pval2)}, ...
        'location', 'southeast');
        l.Position(2) = l.Position(2) - 0.15;
        legend boxoff;
        title(datasetnames{d}); xlabel(parameters{pa});
        vline(0);
        offsetAxes_y;
      end
      if cnt == 4, cnt = 5 ; end
        cnt = cnt + 4;
      end
      print(gcf, '-depsc', sprintf('~/Data/serialHDDM/figure2_posteriors_%s_regression.eps', plots{p}));
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
