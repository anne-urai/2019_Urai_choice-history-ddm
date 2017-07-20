function e8_serialBias_SfN_PPC
% plot posterior predictive checks

addpath(genpath('~/code/Tools'));
warning off; close all; clear;
global datasets datasetnames

set(groot, 'defaultaxesfontsize', 5, 'defaultaxestitlefontsizemultiplier', 1, ...
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

    close all;
    cnt = 1;
    for d = 1:length(datasets),

      try
        % get traces for the model with pupil and rt modulation
        ppc = readtable(sprintf('~/Data/HDDM/%s/stimcoding_nohist/ppq_data.csv', datasets{d}));
        % ppc = readtable('~/Data/ppq_data.csv');
        % ppc = ppc(:, {'rt', 'rt_sampled'}); % save some memory

        % plot the pupil and RT traces
        blues = cbrewer('seq', 'Blues', 3);

        subplot(4,4,cnt); hold on; cnt = cnt + 1;
        h1 = histogram_smooth(ppc.rt, ppc.rt_sampled, blues(1, :), [0 0 0]);

        axis tight; axis square;
        title(datasetnames{d}); xlabel('RT (s)');
        offsetAxes_y;
      end
    end
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/PPC_%s.pdf', plots{p}));
end
end

function h = histogram_smooth(x1, x2, color1, color2)
[f,xi] = ksdensity(x1);
area(xi, f, 'edgecolor', 'none', 'facecolor', color1);
% then the line
[f,xi] = ksdensity(x2);
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
