function barplots_DIC_regression()

% Code to fit the history-dependent drift diffusion models described in
% Urai AE, Gee JW de, Donner TH (2018) Choice history biases subsequent evidence accumulation. bioRxiv:251595
%
% MIT License
% Copyright (c) Anne Urai, 2018
% anne.urai@gmail.com

addpath(genpath('~/code/Tools'));
warning off; close all;
global datasets datasetnames mypath colors

colors(3, :) = mean(colors([1 2], :));
colors = repmat(colors, 7, 1);

% ============================================ %
% DIC COMPARISON BETWEEN DC, Z AND BOTH
% ============================================ %

% 1. STIMCODING, only prevresps
mdls = {'regress_nohist', ...
    'regress_z_lag1', ...
    'regress_dc_lag1', ...
    'regress_dcz_lag1', ...
    'regress_z_lag2', ...
    'regress_dc_lag2', ...
    'regress_dcz_lag2', ...
    'regress_z_lag3', ...
    'regress_dc_lag3', ...
    'regress_dcz_lag3', ...
    'regress_z_lag4', ...
    'regress_dc_lag4', ...
    'regress_dcz_lag4', ...
    'regress_z_lag5', ...
    'regress_dc_lag5', ...
    'regress_dcz_lag5', ...
    'regress_z_lag6', ...
    'regress_dc_lag6', ...
    'regress_dcz_lag6', ...
    'regress_z_lag7', ...
    'regress_dc_lag7', ...
    'regress_dcz_lag7'};

for d = 1:length(datasets),
    close all;
    subplot(4, 4, 1);
    getPlotDIC(mdls, d, colors);
    title(datasetnames{d});
    set(gca, 'xtick', 2:3:length(mdls), 'xticklabel', 1:7);
    ylabel({'\Delta DIC from model'; 'without history'}, 'interpreter', 'tex');
    xlabel('Lag (# trials)')
    drawnow; tightfig;
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/figure1b_HDDM_DIC_regression_d%d.pdf', d));
    fprintf('~/Data/serialHDDM/figure1b_HDDM_DIC_regression_d%d.pdf \n', d)
end

close all;

end

% ============================================ %
% DIC COMPARISON BETWEEN DC, Z AND BOTH
% ============================================ %

function getPlotDIC(mdls, d, colors)

global datasets mypath

axis square; hold on;

mdldic = nan(1, length(mdls));
for m = 1:length(mdls),
    
    if ~exist(sprintf('%s/summary/%s/%s_all.mat', ...
            mypath, datasets{d}, mdls{m}), 'file'),
        disp('cant find this model')
        continue;
    end
    
    load(sprintf('%s/summary/%s/%s_all.mat', ...
        mypath, datasets{d}, mdls{m}));
    
    if (isnan(dic.full) || isempty(dic.full)) && ~all(isnan(dic.chains)),
        dic.full = nanmean(dic.chains);
    end
    mdldic(m) = dic.full;
end

% everything relative to the full model
mdldic = bsxfun(@minus, mdldic, mdldic(1));
mdldic = mdldic(2:end);
[~, bestMdl] = min(mdldic);

for i = 1:length(mdldic),

    if contains(mdls{i+1}, '_z_'),
        xpos = i+0.15;
    elseif contains(mdls{i+1}, '_dc_'),
        xpos = i;
    elseif contains(mdls{i+1}, '_dcz_'), 
        xpos = i-0.15;
    end

    b = bar(xpos, mdldic(i), 'facecolor', colors(i, :), 'barwidth', 0.7, 'BaseValue', 0, ...
        'edgecolor', 'none');
end

% %# Add a text string above/below each bin
% for i = 1:length(mdldic),
%     if mdldic(i) < 0,
%         text(i, mdldic(i) + 0.12*range(get(gca, 'ylim')), ...
%             num2str(round(mdldic(i))), ...
%             'VerticalAlignment', 'top', 'FontSize', 4, 'horizontalalignment', 'center', 'color', 'w');
%     elseif mdldic(i) > 0,
%         text(i, mdldic(i) + 0.12*range(get(gca, 'ylim')), ...
%             num2str(round(mdldic(i))), ...
%             'VerticalAlignment', 'top', 'FontSize', 4, 'horizontalalignment', 'center', 'color', 'w');
%     end
% end

axis tight;
xlim([0.5 length(mdldic)+0.5]);
offsetAxes; box off;
set(gca, 'color', 'none');
set(gca, 'xcolor', 'k', 'ycolor', 'k');

end
