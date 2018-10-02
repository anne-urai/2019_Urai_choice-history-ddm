function barplots_BIC()

% Code to fit the history-dependent drift diffusion models described in
% Urai AE, Gee JW de, Donner TH (2018) Choice history biases subsequent evidence accumulation. bioRxiv:251595
%
% MIT License
% Copyright (c) Anne Urai, 2018
% anne.urai@gmail.com

addpath(genpath('~/code/Tools'));
warning off; close all;
global datasets datasetnames mypath colors

types = {'stimcoding'};
for s = 1:length(types),
    
    % ============================================ %
    % DIC COMPARISON BETWEEN DC, Z AND BOTH
    % ============================================ %
    
    % 1. STIMCODING, only prevresps
    mdls = {'z_prevresp', 'dc_prevresp', ...
        'dc_z_prevresp', 'nohist'};
    for d = 1:length(datasets),
        close all;
        subplot(4, 6, 1);
        getPlotDIC(mdls, d);
         title(datasetnames{d});
        set(gca, 'xtick', 1:3, 'xticklabel', {'z_{bias}', 'v_{bias}', 'both'});
        
        % if ismember(d, [1]),
        ylabel({'\DeltaBIC from model'; 'without history'}, 'interpreter', 'tex');
        % end
        drawnow; tightfig;
        print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/BIC_prevresp_d%d.pdf',  d));
    end
end

close all;

end

% ============================================ %
% DIC COMPARISON BETWEEN DC, Z AND BOTH
% ============================================ %

function getPlotDIC(mdls, d)

global datasets mypath colors
colors(3, :) = mean(colors([1 2], :));
axis square; hold on;

allres = readtable(sprintf('%s/summary/%s/allindividualresults_Gsq.csv', ...
    mypath, datasets{d}));
allres = allres(allres.session == 0, :);
mdldic = nan(size(allres, 1), length(mdls));

for m = 1:length(mdls),
    try
        mdldic(:, m) = allres.(['bic__stimcoding' regexprep(mdls{m}, '_', '')]);
    catch
        continue;
    end
end
assert(~any(isnan(mdldic(:,end))), 'no BICs found');

% everything relative to the nohist model
mdldic = bsxfun(@minus, mdldic, mdldic(:, end));
mdldic = mdldic(:, 1:end-1);

% sum over observers
mdldic = nanmean(mdldic);

% nice looking bargraphs
% colors = [141 165 8;  8 141 165; 150 150 150] ./ 256;
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
            'VerticalAlignment', 'top', 'FontSize', 4, 'horizontalalignment', 'center');
    elseif mdldic(i) > 0,
        text(i, mdldic(i) - 0.12*range(get(gca, 'ylim')), ...
            num2str(round(mdldic(i))), ...
            'VerticalAlignment', 'top', 'FontSize', 4, 'horizontalalignment', 'center', 'color', 'w');
    end
end
xlim([0.5 length(mdldic)+0.5]);
offsetAxes; box off;
axis square;
set(gca, 'color', 'none');

% add text to say how many people win this model?

end
