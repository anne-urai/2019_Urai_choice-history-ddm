
function previousResponseOutcome_DIC

addpath(genpath('~/code/Tools'));
warning off; close all;
global datasets datasetnames mypath

    % ============================================ %
    % DIC COMPARISON BETWEEN DC, Z AND BOTH
    % ============================================ %
    
    % 1. STIMCODING, only prevresps
    mdls = {'dc_z_prevresp', ...
        'dc_prevcorrect', 'nohist'};
    for d = 1:length(datasets),
        close all;
        subplot(4, 5, 1);
        getPlotDIC(mdls, 'stimcoding', d);
        title(datasetnames{d});
        set(gca, 'xtick', 1:2, 'xticklabel', {'dc:c z:c', 'dc:c dc:o'}, 'xticklabelrotation', -30);
        
        %if ismember(d, [1 4]),
            ylabel({'\Delta DIC from model'; 'without history'}, 'interpreter', 'tex');
        %else
        %   ylabel({' '; ' '});
        %end
        drawnow; tightfig;
        print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/DIC_previousResponseOutcome_d%d.pdf', d));
    end

close all;

end

% ============================================ %
% DIC COMPARISON BETWEEN DC, Z AND BOTH
% ============================================ %

function getPlotDIC(mdls, s, d)

global datasets mypath colors
colors(3, :) = mean(colors([1 2], :));
usecolors = [colors(3, :) ; [ 0.984313725490196         0.705882352941177         0.682352941176471]];

axis square; hold on;

mdldic = nan(1, length(mdls));
for m = 1:length(mdls),
    if ~exist(sprintf('%s/summary/%s/%s_%s_all.mat', ...
            mypath, datasets{d}, s, mdls{m}), 'file'),
        disp('cant find this model')
        assert(1==0)
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
    b = bar(i, mdldic(i), 'facecolor', usecolors(i, :), 'barwidth', 0.6, 'BaseValue', 0, ...
        'edgecolor', 'none');
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
            'VerticalAlignment', 'top', 'FontSize', 4, 'horizontalalignment', 'center');
    end
end
axis square; axis tight; 
xlim([0.5 length(mdldic)+0.5]);
offsetAxes; box off;
set(gca, 'color', 'none');
set(gca, 'xcolor', 'k', 'ycolor', 'k');

end
