function e1_serialBias_SfN_DIC()

addpath(genpath('~/code/Tools'));
warning off; close all;
global datasets datasetnames mypath

types = {'stimcoding'};
for s = 1:length(types),

   % ============================================ %
   % DIC COMPARISON BETWEEN DC, Z AND BOTH
    % ============================================ %

    % 1. STIMCODING, only prevresps
    mdls = {'dc_prevcorrect', 'z_prevcorrect', ...
        'dc_z_prevcorrect', 'nohist'};
    for d = 1:length(datasets),
        close all;
        subplot(4, 4, 1);
        try
        getPlotDIC(mdls, types{s}, d, 1);
        title(datasetnames{d});
        set(gca, 'xtick', 1:3, 'xticklabel', {'dc', 'z', 'both'});

        drawnow; tightfig;
        print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/figure1b_HDDM_DIC_%s_prevcorrect_d%d.pdf', types{s}, d));
        fprintf('~/Data/serialHDDM/figure1b_HDDM_DIC_%s_prevcorrect_d%d.pdf \n', types{s}, d);
      end
    end

    % ============================================ %
    % DIC COMPARISON BETWEEN DC, Z AND BOTH
    % ============================================ %

    % build up models
    % 'regress_nohist' % should go last, will be baseline
    % 'regress_dc_prevresp'
    % 'regress_dc_prevresp_prevstim'
    % regress_dc_prevresp_prevstim_vasessions # add learning effects
    % then add modulation
    % then also add multiple responses into the past

    close all; nrsubpl = length(datasets);
    mdls = {'dc_z_prevresp', ...
        'dc_z_prevresp_prevstim',  ...
        'dc_z_prevcorrect', ...
        'dc_z_prev2correct', ...
        'dc_z_prev3correct',  ...
        'nohist'};

    for d = 1:length(datasets),
        subplot(nrsubpl, nrsubpl, d);
        getPlotDIC(mdls, types{s}, d, 1);
        set(gca, 'xtick', 1:length(mdls)-1);
        title(datasetnames{d});
    end

    switch types{s}
        case 'regress'

            subplot(nrsubpl, nrsubpl, nrsubpl+1);
            text(0, -0.2, {'Regression models', ...
                '[1] v ~ 1 + stimulus + prevresp', ...
                '[2] v ~ 1 + stimulus + prevresp + prevstim', ...
                '[3] v ~ 1 + stimulus + prevresp*prevrt + prevstim*prevrt', ...
                '[4] v ~ 1 + stimulus*session + prevresp*prevrt + prevstim*prevrt + prevresp*prevpupil + prevstim*prevpupil', ...
                '[5] v ~ 1 + stimulus*session + prevresp + prevstim', ...
                '    a ~ 1 + session', ...
                '[6] v ~ 1 + session*stimulus + prevresp*prevrt + prevstim*prevrt', ...
                '     a ~ 1 + session', ...
                '[7] v ~ 1 + session*stimulus + prevresp*prevrt + prevstim*prevrt + prevresp*prevpupil + prevstim*prevpupil', ...
                '     a ~ 1 + session', ...
                }, 'fontsize', 6); axis off;
    end

    tightfig;
    % print(gcf, '-depsc', sprintf('~/Data/serialHDDM/suppfigure1b_HDDM_DIC_allmodels_%s_%s.eps', plots{p}, types{s}));
    figure('color', 'none')
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/suppfigure1b_HDDM_DIC_allmodels_%s.pdf',types{s}));
end

close all;

end

% ============================================ %
% DIC COMPARISON BETWEEN DC, Z AND BOTH
% ============================================ %

function getPlotDIC(mdls, s, d, plotBest)

global datasets
axis square; hold on;

mdldic = nan(1, length(mdls));
for m = 1:length(mdls),
    if ~exist(sprintf('/nfs/aeurai/HDDM/summary/%s/%s_%s_all.mat', ...
            datasets{d}, s, mdls{m}), 'file'),
        disp('cant find this model')
        continue;
    end

    load(sprintf('/nfs/aeurai/HDDM/summary/%s/%s_%s_all.mat', ...
        datasets{d}, s, mdls{m}));

    if (isnan(dic.full) || isempty(dic.full)) && ~all(isnan(dic.chains)),
        dic.full = nanmean(dic.chains);
    end
    mdldic(m) = dic.full;
end

if isnan(mdldic(end)), assert(1==0); end

% everything relative to the full mdoel
mdldic = bsxfun(@minus, mdldic, mdldic(end));
mdldic = mdldic(1:end-1);

colors = linspecer(5);
colors = [colors([4 5], :); [0.7 0.7 0.7]; [0.7 0.7 0.7]; [0.7 0.7 0.7];[0.7 0.7 0.7]];

for i = 1:length(mdldic),
    b = bar(i, mdldic(i), 'facecolor', colors(i, :), 'barwidth', 0.6, 'BaseValue', 0);
end

%# Add a text string above/below each bin
for i = 1:length(mdldic),
    if mdldic(i) < 0,
        text(i, mdldic(i) + 0.12*range(get(gca, 'ylim')), ...
            num2str(round(mdldic(i))), ...
            'VerticalAlignment', 'top', 'FontSize', 4, 'horizontalalignment', 'center');
    elseif mdldic(i) > 0,
        text(i, mdldic(i) + 0.12*range(get(gca, 'ylim')), ...
            num2str(round(mdldic(i))), ...
            'VerticalAlignment', 'top', 'FontSize', 4, 'horizontalalignment', 'center');
    end
end
xlim([0.5 length(mdldic)+0.5]);
ylabel('\Delta DIC (from nohist)', 'interpreter', 'tex');
offsetAxes; box off;
axis square;
set(gca, 'color', 'none');

end
