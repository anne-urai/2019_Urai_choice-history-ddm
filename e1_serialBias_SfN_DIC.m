function e1_serialBias_SfN_DIC()

% Code to fit the history-dependent drift diffusion models described in
% Urai AE, Gee JW de, Donner TH (2018) Choice history biases subsequent evidence accumulation. bioRxiv:251595
%
% MIT License
% Copyright (c) Anne Urai, 2018
% anne.urai@gmail.com

addpath(genpath('~/code/Tools'));
warning off; close all;
global datasets datasetnames mypath

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
        getPlotDIC(mdls, types{s}, d);
        title(datasetnames{d}{1});
        set(gca, 'xtick', 1:3, 'xticklabel', {'z_{bias}', 'v_{bias}', 'Both'});
        
		if ismember(d, [1 4]),
        	ylabel({'\Delta DIC from model'; 'without history'}, 'interpreter', 'tex');
		else
			ylabel({' '; ' '});
		end
        drawnow; tightfig;
        print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/figure1b_HDDM_DIC_%s_prevresp_d%d.pdf', types{s}, d));
        fprintf('~/Data/serialHDDM/figure1b_HDDM_DIC_%s_prevresp_d%d.pdf \n', types{s}, d);
    end
    
end

close all;

end

% ============================================ %
% DIC COMPARISON BETWEEN DC, Z AND BOTH
% ============================================ %

function getPlotDIC(mdls, s, d)

global datasets mypath colors
colors(3, :) = mean(colors([1 2], :));
axis square; hold on;

mdldic = nan(1, length(mdls));
for m = 1:length(mdls),
    if ~exist(sprintf('%s/summary/%s/%s_%s_all.mat', ...
            mypath, datasets{d}, s, mdls{m}), 'file'),
        disp('cant find this model')
        continue;
    end
    
    load(sprintf('%s/summary/%s/%s_%s_all.mat', ...
        mypath, datasets{d}, s, mdls{m}));
    
    if (isnan(dic.full) || isempty(dic.full)) && ~all(isnan(dic.chains)),
        dic.full = nanmean(dic.chains);
    end
    mdldic(m) = dic.full;
end

% everything relative to the full model
mdldic = bsxfun(@minus, mdldic, mdldic(end));
mdldic = mdldic(1:end-1);
[~, bestMdl] = min(mdldic);

for i = 1:length(mdldic),
    b = bar(i, mdldic(i), 'facecolor', colors(i, :), 'barwidth', 0.6, 'BaseValue', 0, ...
        'edgecolor', 'none');
end

% [ptchs,ptchGrp] = createPatches(i+1,mdldic(end),0.3, colors(2, :), 0, 0.5);
% hatch(ptchs, [0 3 1], colors(1, :));
% ptchs.EdgeColor = [0 0 0];
% fill the last one

%# Add a text string above/below each bin
for i = 1:length(mdldic),
    if mdldic(i) < 0,
        text(i, mdldic(i) + 0.12*range(get(gca, 'ylim')), ...
            num2str(round(mdldic(i))), ...
            'VerticalAlignment', 'top', 'FontSize', 4, 'horizontalalignment', 'center', 'color', 'w');
    elseif mdldic(i) > 0,
        text(i, mdldic(i) + 0.12*range(get(gca, 'ylim')), ...
            num2str(round(mdldic(i))), ...
            'VerticalAlignment', 'top', 'FontSize', 4, 'horizontalalignment', 'center');
    end
end
axis square; axis tight; 
xlim([0.5 length(mdldic)+0.5]);
offsetAxes; box off;
set(gca, 'color', 'none');
set(gca, 'xcolor', 'k', 'ycolor', 'k');s

end
