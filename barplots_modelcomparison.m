function barplots_modelcomparison()

% Code to fit the history-dependent drift diffusion models described in
% Urai AE, Gee JW de, Donner TH (2018) Choice history biases subsequent evidence accumulation. bioRxiv:251595
%
% MIT License
% Copyright (c) Anne Urai, 2018
% anne.urai@gmail.com

addpath(genpath('~/code/Tools'));
warning off; close all;
global datasets datasetnames mypath

modelIC = {'aic'};
for s = 1:length(modelIC),
    
    % ============================================ %
    % DIC COMPARISON BETWEEN DC, Z AND BOTH
    % ============================================ %
    
    % 1. STIMCODING, only prevresps
    mdls = {'z_prevresp', 'dc_prevresp', ...
        'dc_z_prevresp', 'nohist'};
    for d = 1:length(datasets),
        close all;
        subplot(4, 5, 1);
        getPlotModelIC(mdls, modelIC{s}, d);
        title(datasetnames{d});
        set(gca, 'xtick', 1:3, 'xticklabel', {'z', 'v_{bias}', 'both'});
        ylabel({'\DeltaAIC from model'; 'without history'}, 'interpreter', 'tex');

        drawnow; tightfig;
        print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/modelcomparison_%s_d%d.pdf', modelIC{s}, d));
        disp(d);
    end
end

close all;

end

% ============================================ %
% DIC COMPARISON BETWEEN DC, Z AND BOTH
% ============================================ %

function getPlotModelIC(mdls, s, d)

global datasets mypath colors
colors(3, :) = mean(colors([1 2], :));
axis square; hold on;

mdldic = nan(1, length(mdls));
for m = 1:length(mdls),
    try
    modelcomp = readtable(sprintf('%s/%s/stimcoding_%s/model_comparison.csv', ...
        mypath, datasets{d}, mdls{m}));
    mdldic(m) = modelcomp.(s);
catch
    fprintf('%s/%s/stimcoding_%s/model_comparison.csv  NOT FOUND\n', ...
        mypath, datasets{d}, mdls{m})
    end
end

% everything relative to the full model
mdldic = bsxfun(@minus, mdldic, mdldic(end));
mdldic = mdldic(1:end-1);
[~, bestMdl] = min(mdldic);

for i = 1:length(mdldic),
    if i == bestMdl,
    b = bar(i, mdldic(i), 'facecolor', colors(i, :), 'barwidth', 0.6, 'BaseValue', 0, ...
        'edgecolor', 'k');
    else
    b = bar(i, mdldic(i), 'facecolor', colors(i, :), 'barwidth', 0.6, 'BaseValue', 0, ...
        'edgecolor', 'none');
    end
end

%# Add a text string above/below each bin
for i = 1:length(mdldic),
    if mdldic(i) < 0,
        text(i, mdldic(i) + 0.12*range(get(gca, 'ylim')), ...
            num2str(round(mdldic(i))), ...
            'VerticalAlignment', 'top', 'FontSize', 4, 'horizontalalignment', 'center', 'color', 'w');
    elseif mdldic(i) > 0,
        text(i, mdldic(i) + 0.12*range(get(gca, 'ylim')), ...
            num2str(round(mdldic(i))), ...
            'VerticalAlignment', 'top', 'FontSize', 4, 'horizontalalignment', 'center', 'color', 'w');
    end
end
axis square; axis tight; 
xlim([0.5 length(mdldic)+0.5]);
offsetAxes; box off;
set(gca, 'color', 'none');
set(gca, 'xcolor', 'k', 'ycolor', 'k');

end
