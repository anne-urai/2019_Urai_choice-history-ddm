function e1_serialBias_SfN_DIC()

addpath(genpath('~/code/Tools'));
warning off; close all;
global datasets datasetnames
datasets = {'RT_RDK', 'projects/0/neurodec/Data/MEG-PL', 'NatComm', 'Anke_2afc_neutral', 'Anke_2afc_sequential'};
datasetnames = {'RT', '2IFC',  'NatComm', 'Anke neutral', 'Anke all'};

set(groot, 'defaultaxesfontsize', 7, 'defaultaxestitlefontsizemultiplier', 1.1, ...
    'defaultaxestitlefontweight', 'bold', ...
    'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters');
nrsubpl = length(datasets) + 1;

% ============================================ %
% DIC COMPARISON BETWEEN DC, Z AND BOTH
% ============================================ %

types = {'stimcoding', 'regress'};

for s = 1:2,
    % 1. STIMCODING, only prevresp
    clf;
    mdls = {'dc_prevresp_prevstim', 'z_prevresp_prevstim', ...
        'dc_z_prevresp_prevstim', 'nohist'};
    for d = 1:3,
        subplot(4, 4, d);
        getPlotDIC(mdls, types{s}, d, 1);
        title(['Data: ' datasetnames{d}]);
        set(gca, 'xtick', 1:3, 'xticklabel', {'dc', 'z', 'both'});
    end
    
    subplot(4,4,4);
    switch types{s}
        case 'stimcoding'
            text(-0.2, 0.5, {'Stimcoding models' 'prevresp'}, 'fontsize', 6);
        case 'regress'
            text(-0.2, 0.6, {'Regression models', ...
                'v ~ 1 + stimulus + prevresp', ...
                'z ~ 1 + prevresp'}, 'fontsize', 6);
    end
    axis off;
    print(gcf, '-depsc', sprintf('~/Data/serialHDDM/fig1_DIC_%s.eps', types{s}));
    
end

% % ============================================ %
% DIC COMPARISON BETWEEN DC, Z AND BOTH
% ============================================ %

clf; nrsubpl = 4;
mdls = {'dc_prevresp', 'dc_prevresp_prevstim',  ...
    'dc_prevresp_prevstim_sessions', 'dc_prevresp_prevstim_prevrt', ...
    'dc_prevresp_prevstim_prevrt_prevpupil', ...
    'dc_prevresp_prevstim_prevrt_prevpupil_sessions', 'nohist'};

for d = 1:length(datasets)-2,
    subplot(nrsubpl, nrsubpl, d);
    getPlotDIC(mdls, types{s}, d, 1);
    set(gca, 'xtick', 1:6, 'xticklabel',...
        {'[1]', '[2]', '[3]', '[4]', '[5]', '[6]'});
    title(['Data: ' datasetnames{d}]);
end
subplot(nrsubpl, nrsubpl, 5);
text(0, 0.2, {'Regression models', ...
    '[1] v ~ 1 + stimulus + prevresp', ...
    '[2] v ~ 1 + stimulus + prevresp + prevstim', ...
    '[3] v ~ 1 + stimulus*session + prevresp + prevstim', ...
    '     a ~ 1 + session', ...
    '[4] v ~ 1 + stimulus*session + prevresp*prevrt + prevstim*prevrt', ...
    '     a ~ 1 + session', ...
    '[5] v ~ 1 + stimulus*session + prevresp*prevrt + prevstim*prevrt + prevresp*prevpupil + prevstim*prevpupil', ...
    '     a ~ 1 + session', ...
    '[6] v ~ 1 + session*(stimulus + prevresp*prevrt + prevstim*prevrt + prevresp*prevpupil + prevstim*prevpupil)', ...
    '     a ~ 1 + session', ...
    }, 'fontsize', 7); axis off;
axis off;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/fig1_DIC_prevresp_prevstim.pdf'));

% % ============================================ %
% DIC COMPARISON BETWEEN DC, Z AND BOTH
% ============================================ %

clf;
datasets = {'Anke_2afc_neutral', 'Anke_2afc_repetitive', 'Anke_2afc_alternating'};
datasetnames = {'Anke neutral', 'Anke repetitive', 'Anke alternating'};

mdls = {'dc_prevresp_prevstim', 'z_prevresp_prevstim', ...
    'dc_z_prevresp_prevstim', 'nohist'};

for d = 1:length(datasets),
    subplot(4, 4, d);
    getPlotDIC(mdls, types{s}, d, 1);
    set(gca, 'xtick', 1:3, 'xticklabel',...
        {'dc', 'z', 'both'});
    title(['Data: ' datasetnames{d}]);
end
print(gcf, '-depsc', sprintf('~/Data/serialHDDM/fig1_DIC_biased.eps'));

end

% ============================================ %
% DIC COMPARISON BETWEEN DC, Z AND BOTH
% ============================================ %

function getPlotDIC(mdls, s, d, plotBest)
global datasets
axis square;

mdldic = nan(1, length(mdls));
for m = 1:length(mdls),
    if ~exist(sprintf('~/Data/%s/HDDM/summary/%s_%s_all.mat', ...
            datasets{d}, s, mdls{m}), 'file'),
        continue;
    end
    
    load(sprintf('~/Data/%s/HDDM/summary/%s_%s_all.mat', ...
        datasets{d}, s, mdls{m}));
    
    if (isempty(dic.full) || isnan(dic.full)) && ~all(isnan(dic.chains)),
        dic.full = nanmean(dic.chains);
    end
    mdldic(m) = dic.full;
    allfull(m, :) = dic.chains;
end

% if there are two, take the mean
if numel(mdldic) == 3 || isnan(mdldic(4)),
    mdldic(4) = nanmean(mdldic(1:3));
    mdldic(4) = 290618;
end

% everything relative to the last model
mdldic = bsxfun(@minus, mdldic, mdldic(end));
mdldic = mdldic(1:end-1);

% now plot
bar(mdldic, 'facecolor', linspecer(1), 'barwidth', 0.4);

%# Add a text string above/below each bin
for i = 1:length(mdldic),
    if mdldic(i) < 0,
        text(i, mdldic(i) - 0.06*range(get(gca, 'ylim')), ...
            num2str(round(mdldic(i))), ...
            'VerticalAlignment', 'top', 'FontSize', 4, 'horizontalalignment', 'center');
    elseif mdldic(i) > 0,
        text(i, mdldic(i) + 0.12*range(get(gca, 'ylim')), ...
            num2str(round(mdldic(i))), ...
            'VerticalAlignment', 'top', 'FontSize', 4, 'horizontalalignment', 'center');
    end
end

% indicate which one is best
bestcolor = linspecer(3);
if plotBest,
    [~, idx] = min(mdldic);
    hold on;
    disp(idx);
    bar(idx, mdldic(idx), 'facecolor', bestcolor(3, :), 'barwidth', 0.4);
end

box off;
ylabel('\Delta DIC (from nohist)');
set(gca, 'xticklabel', {'dc', 'z'});
% set(gca, 'YTickLabel', num2str(get(gca, 'YTick')'));
xlim([0 length(mdls)-0.5]);
axis square;
ylim(1.1*get(gca, 'ylim'));

end