function e1_serialBias_SfN_DIC()

addpath(genpath('~/code/Tools'));
warning off; close all;
global datasets datasetnames
datasets = {'RT_RDK', 'projects/0/neurodec/Data/MEG-PL', 'NatComm', 'Anke_2afc_neutral', 'Anke_2afc_sequential'};
datasetnames = {'RT', '2IFC',  'NatComm', 'Anke neutral', 'Anke all'};

set(groot, 'defaultaxesfontsize', 6, 'defaultaxestitlefontsizemultiplier', 1.1, ...
    'defaultaxestitlefontweight', 'bold', ...
    'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters');
nrsubpl = length(datasets) + 1;

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
        subplot(nrsubpl, nrsubpl,d);
        getPlotDIC(mdls, types{s}, d);
        title(['Data: ' datasetnames{d}]);
    end
    
    subplot(nrsubpl,nrsubpl,d+1);
    switch types{s}
        case 'stimcoding'
            text(-0.2, 0.5, {'Stimcoding models' 'prevresp'}, 'fontsize', 6);
        case 'regress'
            text(-0.2, 0.6, {'Regression models', ...
                'v ~ 1 + stimulus + prevresp', ...
                'z ~ 1 + prevresp'}, 'fontsize', 6);
    end
    axis off;
    
    % 2. STIMCODING, PREVRESP + PREVSTIM
    mdls = {'dc_prevresp_prevstim', 'z_prevresp_prevstim', ...
        'dc_z_prevresp_prevstim'};
    for d = 1:length(datasets),
        subplot(nrsubpl,nrsubpl,d+nrsubpl);
        getPlotDIC(mdls, types{s}, d);
    end
    subplot(nrsubpl,nrsubpl,d+nrsubpl+1);
    switch types{s}
        case 'stimcoding'
            text(-0.2, 0.5, {'Stimcoding models', 'prevresp + prevstim'}, 'fontsize', 6);
        case 'regress'
            text(-0.2, 0.6, {'Regression models', ...
                'v ~ 1 + stimulus + prevresp + prevstim', ...
                'z ~ 1 + prevresp + prevstim'}, 'fontsize', 6);
    end
    axis off;
    print(gcf, '-depsc', sprintf('~/Data/serialHDDM/fig1_DIC_%s.eps', types{s}));
    
end

% ============================================ %
% DIC COMPARISON BETWEEN DC, Z AND BOTH
% ============================================ %

clf; nrsubpl = length(datasets);
mdls = {'dc_prevresp', 'dc_prevresp_prevstim', 'dc_prevresp_prevstim_prevrt', ...
    'dc_prevresp_prevstim_prevrt_prevpupil', 'dc_prevresp_prevstim_prevrt_prevpupil_sessions', ...
    'dc_prevresp_prevstim_sessions'};
for d = 1:length(datasets),
    subplot(nrsubpl, nrsubpl, d);
    getPlotDIC(mdls, types{s}, d);
    set(gca, 'xtick', 1:5, 'xticklabel',...
        {'[1]', '[2]', '[4]', '[5]', '[6]'});
    title(['Data: ' datasetnames{d}]);
    ylabel('\Delta DIC (from [3])');
end
subplot(nrsubpl,nrsubpl,d+3);
text(-0.2, 0.6, {'Regression models', ...
    '[1] v ~ 1 + stimulus + prevresp', ...
    '[2] ... + prevstim', ...
    '[3] ... + v:session + a:session', ...
    '[4] ... + prevrt*prevresp + prevrt*prevstim', ...
    '[5] ... + prevpupil*prevresp + prevpupil*prevstim', ...
    '[6] ... + serialbias:session', ...
    }, 'fontsize', 6); axis off;
subplot(5,5,11); plot(1,1, 'w.'); axis off;
print(gcf, '-depsc', sprintf('~/Data/serialHDDM/fig1_DIC_prevresp_prevstim.eps'));

end

% ============================================ %
% DIC COMPARISON BETWEEN DC, Z AND BOTH
% ============================================ %

function getPlotDIC(mdls, s, d)
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
    
    if isnan(dic.full) && ~all(isnan(dic.chains)),
        dic.full = nanmean(dic.chains);
    end
    mdldic(m) = dic.full;
end

% if there are two, take the mean
if numel(mdldic) == 2 | isnan(mdldic(3)), 
    mdldic(3) = nanmean(mdldic(1:2));
end

% everything relative to the full mdoel
mdldic = bsxfun(@minus, mdldic, mdldic(end));
bar(mdldic(1:end-1), 'facecolor', linspecer(1));

%# Add a text string above/below each bin
for i = 1:length(mdldic)-1,
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

box off;
ylabel('\Delta DIC (from full)');
set(gca, 'xticklabel', {'dc', 'z'});
% set(gca, 'YTickLabel', num2str(get(gca, 'YTick')'));
xlim([0 length(mdls)]);
axis square;

end