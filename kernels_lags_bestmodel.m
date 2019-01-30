function kernels_lags_bestmodel

global mypath datasets 
addpath(genpath('~/code/Tools'));
warning off; close all;

% ALL MODELS THAT WERE RAN
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
    'regress_dc_lag7-10', ...
    'regress_z_lag7-10',...
    'regress_dcz_lag7-10',...
    'regress_dc_lag11-15', ...
    'regress_z_lag11-15',...
    'regress_dcz_lag11-15'};

numlags = 8;
lagnames = {'1', '2', '3', '4', '5', '6', '7-10', '11-15'};
vars = {'z_correct', 'z_error', 'v_correct', 'v_error'};
for m = 1:length(vars),
    alldata.(vars{m})   = nan(length(datasets), numlags);
    alldata.([vars{m} '_pval'])   = nan(length(datasets), numlags);
end

mat_z.r     = nan(length(datasets), numlags);
mat_z.pval  = nan(length(datasets), numlags);
mat_dc.r    = nan(length(datasets), numlags);
mat_dc.pval = nan(length(datasets), numlags);

for d = 1:length(datasets),
    
    % ============================= %
    % 1. DETERMINE THE BEST MODEL
    % ============================= %
    
    mdldic = nan(1, length(mdls));
    for m = 1:length(mdls),
        
        if ~exist(sprintf('%s/summary/%s/%s_all.mat', ...
                mypath, datasets{d}, mdls{m}), 'file'),
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
    bestmodelname = regexprep(mdls{bestMdl+1}, '_', '');
    
    % ========================================================== %
    % use the full model with both starting point and drift bias for plots
    % ========================================================== %

    useFullModel = false;
    if useFullModel,
        bestMdl = 9;
        bestmodelname = 'regressdczlag3';
    end
    
    % ========================================================== %
    % 2. FOR THIS MODEL, RECODE INTO CORRECT AND ERROR
    % ========================================================== %
    
    dat = readtable(sprintf('%s/summary/%s/allindividualresults.csv', mypath, datasets{d}));
    traces = readtable(sprintf('%s/%s/%s/group_traces.csv', mypath, datasets{d}, mdls{bestMdl+1}));

        for l = 1:numlags,
            if l == 1,
                lname = [];
            else
                lname = l;
            end
                
        for v = 1:length(vars),
            try
                switch vars{v}
                case 'z_correct'
                alldata.(vars{v})(d,l) = ...
                    nanmean(dat.(['z_prev' num2str(lname) 'resp__' bestmodelname]) + ...
                    dat.(['z_prev' num2str(lname) 'stim__' bestmodelname]));
                alldata.([vars{v} '_pval'])(d,l) = posteriorpval(traces.(['z_prev' num2str(lname) 'resp']) + ...
                    traces.(['z_prev' num2str(lname) 'stim']), 0);

                case 'z_error'
                    alldata.z_error(d,l) = ...
                nanmean(dat.(['z_prev' num2str(lname) 'resp__' bestmodelname]) - ...
                dat.(['z_prev' num2str(lname) 'stim__' bestmodelname]));
                                alldata.([vars{v} '_pval'])(d,l) = posteriorpval(traces.(['z_prev' num2str(lname) 'resp']) - ...
                    traces.(['z_prev' num2str(lname) 'stim']), 0);
                case 'v_correct'
                    alldata.v_correct(d,l) = ...
                    nanmean(dat.(['v_prev' num2str(lname) 'resp__' bestmodelname]) + ...
                    dat.(['v_prev' num2str(lname) 'stim__' bestmodelname]));
                                    alldata.([vars{v} '_pval'])(d,l) = posteriorpval(traces.(['v_prev' num2str(lname) 'resp']) + ...
                    traces.(['v_prev' num2str(lname) 'stim']), 0);
                case 'v_error'
                     alldata.v_error(d,l) = ...
                nanmean(dat.(['v_prev' num2str(lname) 'resp__' bestmodelname]) - ...
                dat.(['v_prev' num2str(lname) 'stim__' bestmodelname]));
                                alldata.([vars{v} '_pval'])(d,l) = posteriorpval(traces.(['v_prev' num2str(lname) 'resp']) - ...
                    traces.(['v_prev' num2str(lname) 'stim']), 0);
            end
        end

        end
    
    end
    
    % ========================================================== %
    % ALSO COMPUTE CORRELATIONS
    % TO DO: PARTIAL OUT THE EFFECT OF PREVIOUS REPETITIONS!
    % ========================================================== %

    %covariates = zeros(size(dat.(['repetition' num2str(l)])));
    for l = 1:numlags,

        if l == 1,
            lname = [];
        else
            lname = l;
        end
        
        repeat = dat.(['repetition_corrected' num2str(l)]);
        %repeat = dat.(['repetition' num2str(l)]);

        try
            [mat_z.r(d, l), mat_z.ci(d,l,:), mat_z.pval(d,l)] = ...
                spearmans(dat.(['z_prev' num2str(lname) 'resp__' bestmodelname]), ...
                repeat);
        end

        %     [mat_z.r(d, l), mat_z.pval(d,l)] = ...
        %         partialcorr(dat.(['z_prev' num2str(lname) 'resp__' bestmodelname]), ...
        %         repeat, covariates, 'type', 'spearman', 'rows', 'complete');
        % end
        try
            [mat_dc.r(d, l), mat_dc.ci(d,l,:), mat_dc.pval(d,l)] = ...
                spearmans(dat.(['v_prev' num2str(lname) 'resp__' bestmodelname]), ...
                repeat);
        end

        % [mat_dc.r(d, l), mat_dc.pval(d,l)] = ...
        %         partialcorri(dat.(['v_prev' num2str(lname) 'resp__' bestmodelname]), ...
        %         repeat, covariates, 'type', 'spearman', 'rows', 'complete');
        % end
        %covariates = [covariates repeat];

        % if l == 1,
        %     covariates = covariates(:, 2:end);
        % elseif l == 4,
        %         assert(1==0)
        % end
    end
    
end

% ========================================================== %
% 3. PLOT THE VARIABLES THAT ARE PRESENT FOR THIS BEST MODEL
% ========================================================== %

colors = cbrewer('qual', 'Set2', length(datasets));

% CREATE FIGURE
for pltidx = 1:length(vars),
    
    close all;
    sp1 = subplot(4,4,1); hold on;
    plot([1 numlags], [0 0], 'k', 'linewidth', 0.5);
    
    for d = 1:length(datasets),
        plot(1:numlags, alldata.(vars{pltidx})(d, :), 'color', colors(d, :), 'linewidth', 1);
    
        h = (alldata.([vars{pltidx} '_pval'])(d,:) < 0.05);
        if any(h>0),
            plot(find(h==1), alldata.(vars{pltidx})(d, (h==1)), '.', 'markeredgecolor', colors(d, :), ...
                'markerfacecolor', colors(d,:), 'markersize', 7);
        end
    end

    % average across datasets
    plot(1:numlags, nanmean(alldata.(vars{pltidx})), 'k', 'linewidth', 1);
    [h, pval] = ttest(alldata.(vars{pltidx}));
    if any(h>0),
        plot(find(h==1), nanmean(alldata.(vars{pltidx})(:, (h==1))), ...
            'k.', 'markersize', 10);
    end
    xlabel('Lags (# trials)');
    ylabel(regexprep(regexprep(vars{pltidx}, '_', ' ~ previous '), 'v ', 'v_{bias} '));
    set(gca, 'xtick', 1:numlags, 'xticklabel', lagnames, 'xticklabelrotation', -20, 'xcolor', 'k', 'ycolor', 'k');
    axis tight; offsetAxes;
    
    tightfig;
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/regressionkernels_correcterror_%d.pdf', pltidx));
    fprintf('~/Data/serialHDDM/regressionkernels_correcterror_%d.pdf \n', pltidx)
end


% ========================================================== %
% 4. PLOT CORRELATION KERNELS
% ========================================================== %


% Z CORRELATION KERNELS
close all;
subplot(441); hold on;
plot([1 numlags], [0 0], 'k', 'linewidth', 0.5);
for d = 1:length(datasets),
    plot(1:numlags, mat_z.r(d, :), 'color', colors(d, :), 'linewidth', 1);
    h = (mat_z.r(d,:) < 0.05);
    if any(h>0),
        plot(find(h==1), mat_z.r(d, (h==1)), '.', 'markeredgecolor', colors(d, :), ...
            'markerfacecolor', colors(d,:), 'markersize', 7);
    end
end

plot(1:numlags, nanmean(mat_z.r), 'k', 'linewidth', 1);
[h, pval] = ttest(mat_z.r);
if any(h>0),
    plot(find(h==1), nanmean(mat_z.r(:, (h==1))), ...
        '.k',  'markersize', 10);
end

ylabel({'Correlation, P(repeat) with' 'z ~ previous response'})
set(gca, 'xtick', 1:numlags, 'xticklabel', lagnames, 'xticklabelrotation', -20, 'xcolor', 'k', 'ycolor', 'k');
xlabel('Lags (# trials)');
axis tight;
ylim([-0.5 1]);
offsetAxes; tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/correlationkernels_z.pdf'));

% DC CORRELATION KERNELS
close all;
subplot(441); hold on;
plot([1 numlags], [0 0], 'k', 'linewidth', 0.5);
for d = 1:length(datasets),
    plot(1:numlags, mat_dc.r(d, :), 'color', colors(d, :), 'linewidth', 1);
    h = (mat_dc.pval(d,:) < 0.05);
    if any(h>0),
        plot(find(h==1), mat_dc.r(d, (h==1)), '.', 'markeredgecolor', colors(d, :), ...
            'markerfacecolor', colors(d,:), 'markersize', 7);
    end
end

plot(1:numlags, nanmean(mat_dc.r), 'k', 'linewidth', 1);
[h, pval] = ttest(mat_dc.r);
if any(h>0),
    plot(find(h==1), nanmean(mat_dc.r(:, (h==1))), ...
        '.k',  'markersize', 10);
end

ylabel({'Correlation, P(repeat) with' 'v_{bias} ~ previous response'})
set(gca, 'xtick', 1:numlags, 'xticklabel', lagnames, 'xticklabelrotation', -20, 'xcolor', 'k', 'ycolor', 'k');
xlabel('Lags (# trials)');
axis tight;
ylim([-0.5 1]);
offsetAxes; tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/correlationkernels_dc.pdf'));

end
