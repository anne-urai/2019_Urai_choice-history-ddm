function e1_serialBias_SfN_DIC()

addpath(genpath('~/code/Tools'));
warning off; close all;
global datasets datasetnames
datasets = {'RT_RDK', 'projects/0/neurodec/Data/MEG-PL', 'Anke_2afc_sequential'};
datasetnames = {'RT', '2IFC', 'Anke'};

set(groot, 'defaultaxesfontsize', 7, 'defaultaxestitlefontsizemultiplier', 1.1, ...
    'defaultaxestitlefontweight', 'bold', ...
    'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters');

% ============================================ %
% DIC COMPARISON BETWEEN DC, Z AND BOTH
% ============================================ %

types = {'stimcoding', 'regress'};

for s = 1:2
    % 1. STIMCODING, only prevresp
    clf;
    mdls = {'dc_prevresp', 'z_prevresp', ...
        'dc_z_prevresp'};
    for d = 1:length(datasets),
        subplot(4,4,d);
        getPlotDIC(mdls, types{s}, d);
        title(['Data: ' datasetnames{d}]);
    end
    
    subplot(4,4,d+1);
    switch types{s}
        case 'stimcoding'
            text(0, 0.5, {'Stimcoding models' 'prevresp'}, 'fontsize', 6);
        case 'regress'
            text(0, 0.7, {'Regression models', ...
                'v ~ 1 + stimulus + prevresp', ...
                'z ~ 1 + prevresp'}, 'fontsize', 6);
    end
    axis off;
    
    % 2. STIMCODING, PREVRESP + PREVSTIM
    mdls = {'dc_prevresp_prevstim', 'z_prevresp_prevstim', ...
        'dc_z_prevresp_prevstim'};
    for d = 1:length(datasets),
        subplot(4,4,d+4);
        getPlotDIC(mdls, types{s}, d);
    end
    subplot(4,4,d+5);
    switch types{s}
        case 'stimcoding'
            text(0, 0.5, {'Stimcoding models' 'prevresp + prevstim'}, 'fontsize', 6);
        case 'regress'
            text(0, 0.7, {'Regression models', ...
                'v ~ 1 + stimulus + prevresp + prevstim', ...
                'z ~ 1 + prevresp + prevstim'}, 'fontsize', 6);
    end
    axis off;
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/fig1_DIC_%s.pdf', types{s}));
    
end


end
% ============================================ %
% DIC COMPARISON BETWEEN DC, Z AND BOTH
% ============================================ %

function getPlotDIC(mdls, s, d)
global datasets

mdldic = nan(1, length(mdls));
for m = 1:length(mdls),
    if ~exist(sprintf('~/Data/%s/HDDM/summary/%s_%s_all.mat', ...
            datasets{d}, s, mdls{m}), 'file'),
        return;
    end
    
    load(sprintf('~/Data/%s/HDDM/summary/%s_%s_all.mat', ...
        datasets{d}, s, mdls{m}));
    
    if isnan(dic.full) && ~all(isnan(dic.chains)),
        dic.full = nanmean(dic.chains);
    end
    mdldic(m) = dic.full;
end

% everything relative to the full mdoel
mdldic = bsxfun(@minus, mdldic, mdldic(3));
bar(mdldic(1:2), 'facecolor', linspecer(1));

%# Add a text string above/below each bin
for i = 1:2,
    if mdldic(i) < 0,
        text(i, mdldic(i) - 0.05*range(get(gca, 'ylim')), ...
            num2str(round(mdldic(i))), ...
            'VerticalAlignment', 'top', 'FontSize', 6, 'horizontalalignment', 'center');
    elseif mdldic(i) > 0,
        text(i, mdldic(i) + 0.1*range(get(gca, 'ylim')), ...
            num2str(round(mdldic(i))), ...
            'VerticalAlignment', 'top', 'FontSize', 6, 'horizontalalignment', 'center');
    end
end

box off;
ylabel('\Delta DIC (from full model)');
set(gca, 'xticklabel', {'dc', 'z'});
% set(gca, 'YTickLabel', num2str(get(gca, 'YTick')'));
axis square;
xlim([0 3]);

end