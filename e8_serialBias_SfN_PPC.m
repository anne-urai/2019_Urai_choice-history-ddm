function e8_serialBias_SfN_PPC
% plot posterior predictive checks

addpath(genpath('~/code/Tools'));
warning off; close all; clear;
global datasets datasetnames

% neutral vs biased plots
datasets = {'RT_RDK', 'NatComm', 'MEG', 'Anke_2afc_sequential'};
datasetnames = { {'2AFC, RT', 'n = 22'}, ...
    {'2IFC, Urai et al. 2017', 'n = 27'}, {'2IFC, replication', 'n = 61'}, ...
    {'2AFC, Braun et al. 2017', 'n = 22'}};

% ========================================== %
% MODULATION OF SERIAL CHOICE BIAS
% ========================================== %
for d = [1 2 4 3],
    
    if ~exist(sprintf('~/Data/HDDM/%s/stimcoding_nohist/ppq_data.csv', datasets{d}), 'file'),
        fprintf('cannot find ~/Data/HDDM/%s/stimcoding_nohist/ppq_data.csv \n', datasets{d});
        continue;
    else
        disp(datasets{d});
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
    
    axis tight; axis square; 
    offsetAxes_y;
    xlim([-3 3]); set(gca, 'xtick', [-3 0 3]);
    if d > 1,
        set(gca, 'yticklabel', []);
    else
        ylabel('Probability');
    end
    tightfig;
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/PPC_d%d.pdf', d));
    
    % keep info about the distribution of errors
    errors{d} = ppc.rt - ppc.rt_sampled;
    
end

%% also show a histogram of the rt error
close all;

colors = cbrewer('qual', 'Dark2', length(errors));
subplot(441); hold on;
for e = 1:length(errors),
    h(e) = histogram(errors{e}, 'displaystyle', 'stairs', 'normalization', 'pdf', ...
        'edgecolor', colors(e, :));
end
xlim([-2 2]);
xlabel('$$RT-\widehat{RT}$$','Interpreter','Latex');
ylabel('Probability');

l = legend([h([1 4 2 3])], ...
    {datasetnames{1}{1} datasetnames{4}{1} datasetnames{2}{1} datasetnames{3}{1}}, 'location', 'southeast');
l.Box = 'off';
l.Position(1) = l.Position(1) + 0.45;
l.Position(2) = l.Position(2) - 0.1;

subplot(442); axis off;
tightfig;

print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/PPC_compare.pdf'));


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
h = plot(xi, f, 'color', 'k', 'linewidth', 1);

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
