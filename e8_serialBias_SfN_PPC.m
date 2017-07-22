function e8_serialBias_SfN_PPC
% plot posterior predictive checks

addpath(genpath('~/code/Tools'));
warning off; close all; clear;
global datasets datasetnames

% ========================================== %
% MODULATION OF SERIAL CHOICE BIAS
% ========================================== %
for d = 1:length(datasets),

    if ~exist(sprintf('~/Data/HDDM/%s/stimcoding_nohist/ppq_data.csv', datasets{d}), 'file'),
      fprintf('cannot find ~/Data/HDDM/%s/stimcoding_nohist/ppq_data.csv \n', datasets{d});
        continue;
    end

    close all; subplot(4,4,1); hold on;
    title(datasetnames{d});
    xlabel('RT (s)');

    % get traces for the model with pupil and rt modulation
    ppc = readtable(sprintf('~/Data/HDDM/%s/stimcoding_nohist/ppq_data.csv', datasets{d}));

    % make sure errors are negative
    ppc.correct = (ppc.stimulus == ppc.response);
    ppc.rt(ppc.correct == 1)           = abs(ppc.rt(ppc.correct == 1));
    ppc.rt(ppc.correct == 0)           = -abs(ppc.rt(ppc.correct == 0));
    ppc.rt_sampled(ppc.correct == 1)   = abs(ppc.rt_sampled(ppc.correct == 1));
    ppc.rt_sampled(ppc.correct == 0)   = -abs(ppc.rt_sampled(ppc.correct == 0));
    ppc = ppc(:, {'rt', 'rt_sampled', 'correct'}); % save some memory

    % plot the pupil and RT traces
    bestcolor = linspecer(4, 'qualitative');
    histogram_smooth(ppc.rt, ppc.rt_sampled, bestcolor(3, :), bestcolor(2, :));

    axis tight; axis square; xlim([-4 4]);
    offsetAxes_y; ylabel('Probability');

    %% also show a histogram of the rt error
    sp2 = subplot(445);
    error = (ppc.rt - ppc.rt_sampled);

    % then the line
    violinPlot(error, 'histOpt', 1, 'showMM', 0, 'color', [0.5 0.5 0.5]);
    ylim([-4 4]);
    xlim([0 2]); set(gca, 'xtick', [0.5 1.5], 'xticklabel', []);
    set(gca, 'xcolor', 'w');
    axis square;
    ylabel('$$RT-\widehat{RT}$$','Interpreter','Latex');
    sp2.Position(2) = sp2.Position(2) - 0.01;
    tightfig;

    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/PPC_d%d.pdf', d));
end
end

function h = histogram_smooth(x1, x2, color1, color2)

% manually count so i can plot myself
[n, edges] = histcounts(x1, 100, 'normalization', 'pdf');

posidx = find(edges > 0); posidx(posidx > length(n)) = [];
negidx = find(edges < 0);

bar(edges(posidx), n(posidx), 'edgecolor', 'none', 'facecolor', color1, 'barwidth', 1);
bar(edges(negidx), n(negidx), 'edgecolor', 'none', 'facecolor', color2, 'barwidth', 1);

% then the line
[f,xi] = ksdensity(x2);
h = plot(xi, f, 'color', 'k', 'linewidth', 0.5);
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
