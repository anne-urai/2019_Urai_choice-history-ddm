function e8_serialBias_SfN_PPC
% plot posterior predictive checks

addpath(genpath('~/code/Tools'));
warning off; close all; clear;
global datasets datasetnames mypath

% ========================================== %
% MODULATION OF SERIAL CHOICE BIAS
% ========================================== %

plotWhich = 'error'; % {'error', 'biased'};


for d = 1:length(datasets),
    
    if ~exist(sprintf('%s/%s/stimcoding_nohist/ppc_data.csv', mypath, datasets{d}), 'file'),
        fprintf('cannot find %s/stimcoding_nohist/ppc_data.csv \n', datasets{d});
        continue;
    else
        disp(datasets{d});
    end
    
    close all; subplot(4,4,1); hold on;
    xlabel('RT (s)');
    
    % get traces for the model with pupil and rt modulation
    ppc = readtable(sprintf('%s/%s/stimcoding_nohist/ppc_data.csv', mypath, datasets{d}));
    ppc.correct                        = (ppc.stimulus == ppc.response);
    ppc.repeat                         = zeros(size(ppc.response));
    ppc.repeat(ppc.response == (ppc.prevresp > 0)) = 1;
    
    % for each observers, compute their bias
    [gr, sjs] = findgroups(ppc.subj_idx);
    sjrep = splitapply(@nanmean, ppc.repeat, gr);
    sjrep = sjs(sjrep < 0.5);
    
    % recode real data into biased vs unbiased
    ppc.biased                         = ppc.repeat;
    altIdx                             = ismember(ppc.subj_idx, sjrep);
    ppc.biased(altIdx) = double(~(ppc.biased(altIdx))); % flip
 
    switch plotWhich
        case 'error'
            ppc.biased = ppc.correct;
    end
    
    % unbiased RTs negative
    ppc.rt(ppc.biased == 1)           = abs(ppc.rt(ppc.biased == 1));
    ppc.rt(ppc.biased == 0)           = -abs(ppc.rt(ppc.biased == 0));
    
    % SAME FOR THE SIMULATED DATA
    ppc.correct_sampled = (ppc.stimulus == ppc.response_sampled);
    
    % recode into repeat and alternate for the model
    ppc.repeat_sampled = zeros(size(ppc.response_sampled));
    ppc.repeat_sampled(ppc.response_sampled == (ppc.prevresp > 0)) = 1;
    
    % recode into biased and unbiased choices
    ppc.biased_sampled = ppc.repeat_sampled;
    altIdx = ismember(ppc.subj_idx, sjrep);
    ppc.biased_sampled(altIdx) = double(~(ppc.biased_sampled(altIdx))); % flip
    
    switch plotWhich
        case 'error'
            ppc.biased_sampled = ppc.correct_sampled;
    end
    
    % define the sampled RT also by the sampled correctness!
    ppc.modelcorrect                   = (ppc.response_sampled == ppc.stimulus);
    ppc.rt_sampled(ppc.biased_sampled == 1)   = abs(ppc.rt_sampled(ppc.biased_sampled == 1));
    ppc.rt_sampled(ppc.biased_sampled == 0)   = -abs(ppc.rt_sampled(ppc.biased_sampled == 0));
    ppc = ppc(:, {'rt', 'rt_sampled'}); % save some memory
    
    % plot the pupil and RT traces
    
    switch plotWhich
        case 'error'
            bestcolor = linspecer(4, 'qualitative');
            bestcolor = bestcolor([3 2], :);
        case 'biased'
            bestcolor = cbrewer('div', 'PiYG', 6);
            bestcolor = bestcolor([1 end], :);
    end
    histogram_smooth(ppc.rt, ppc.rt_sampled, bestcolor(1, :), bestcolor(2, :));
    
    axis tight; axis square;
    offsetAxes_y;
	maxRT = round(max(abs(ppc.rt)));
	if maxRT == 5, maxRT = 4; end
    xlim([-maxRT maxRT]); set(gca, 'xtick', [-maxRT 0 maxRT]);
	disp(maxRT);
    ylabel('Probability');
    set(gca, 'yticklabel', []);
    
    tightfig;
    switch plotWhich
        case 'error'
            print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/PPC_d%d.pdf', d));
        case 'biased'
            print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/PPC_d%d_biased.pdf', d));
    end
    
end

%% also show a histogram of the rt error
% 
% ds = [1 6 2 3 9 7];
% close all;
% subplot(441); 
% % https://nl.mathworks.com/matlabcentral/answers/60818-boxplot-with-vectors-of-different-lengths
% col = @(x)reshape(x,numel(x),1);
% boxplot2 = @(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),...
%     cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});
% boxplot2(errors(ds));
% 
% set(gca, 'xtick', 1:length(ds), 'xticklabel', datasetnames{ds});
% ylabel('$$|\widehat{RT}-RT|$$','Interpreter','Latex');
% print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/PPC_compare.pdf'));
% 
% 
% colors = cbrewer('qual', 'Dark2', length(errors));
% subplot(441); hold on;
% for e = 1:length(errors),
%     plot([median(errors{e}) median(errors{e})], [0 4], ':', 'color', colors(e, :), 'linewidth', 0.2);
%     h(e) = histogram(errors{e}, 'displaystyle', 'stairs', 'normalization', 'pdf', ...
%         'edgecolor', colors(e, :));
% end
% xlim([0 2]); ylim([0 4]);
% set(gca, 'yticklabel', []);
% ylabel('Probability');
% tightfig;
% offsetAxes;

% how often does the model make the same choice as the subject?
%disp(estimperf);
%disp(datasetnames);

end

function h = histogram_smooth(x1, x2, color1, color2)

% manually count so i can plot myself
[n, edges] = histcounts(x1, -3:0.1:3, 'normalization', 'pdf');

posidx = find(edges > 0); posidx(posidx > length(n)) = [];
negidx = find(edges < 0);

bar(edges(posidx), n(posidx), 'edgecolor', 'none', 'facecolor', color1, 'barwidth', 1);
bar(edges(negidx), n(negidx), 'edgecolor', 'none', 'facecolor', color2, 'barwidth', 1);

% then the line
[f,xi] = ksdensity(x2);
h = plot(xi, f, 'color', 'k', 'linewidth', 1);
set(gca, 'color', 'none');

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