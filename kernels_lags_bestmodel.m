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
    'regress_dcz_lag6'};

numlags = 6;
lagnames = {'1', '2', '3', '4', '5', '6', '7-10', '11-15'};
vars = {'z_correct', 'z_error', 'v_correct', 'v_error', ...
    'z_prevresp', 'z_prevstim', 'v_prevresp', 'v_prevstim'};
for m = 1:length(vars),
    alldata.(vars{m})       = nan(length(datasets), numlags);
    alldata.([vars{m} '_fullmodel'])       = nan(length(datasets), numlags);
    alldata.([vars{m} '_pval'])   = nan(length(datasets), numlags);
end
fullmodelname = 'regressdczlag6';
global individualrep

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
    bestmodelname = regexprep(regexprep(mdls{bestMdl+1}, '_', ''), '-', 'to');
    bestmodelnames{d} = bestmodelname;
    disp(bestmodelname)

    % ========================================================== %
    % 2. FOR THIS MODEL, RECODE INTO CORRECT AND ERROR
    % ========================================================== %
    
    dat = readtable(sprintf('%s/summary/%s/allindividualresults.csv', mypath, datasets{d}));
    dat = dat(dat.session == 0, :);
    traces = readtable(sprintf('%s/%s/%s/group_traces.csv', mypath, datasets{d}, mdls{bestMdl+1}));
    individualrep = sign(dat.repetition - 0.5);

        for l = 1:numlags,
            if l == 1,
                lname = '';
            elseif l == 7,
                lname = '7_10';
            elseif l == 8,
                lname = '11_15';
            else
                lname = num2str(l);
            end
                
        for v = 1:length(vars),
                switch vars{v}
                case 'z_correct'
                    try
                alldata.([vars{v} '_fullmodel'])(d,l) = ...
                    summarize(dat.(['z_prev' lname 'resp__' fullmodelname]) + ...
                    dat.(['z_prev' lname 'stim__' fullmodelname]));
                    end
                    try
                alldata.(vars{v})(d,l) = ...
                    summarize(dat.(['z_prev' lname 'resp__' bestmodelname]) + ...
                    dat.(['z_prev' lname  'stim__' bestmodelname]));

                alldata.([vars{v} '_pval'])(d,l) = posteriorpval(traces.(['z_prev' lname  'resp']) + ...
                    traces.(['z_prev' lname  'stim']), 0);
                    end

                case 'z_error'
                    try
                alldata.([vars{v} '_fullmodel'])(d,l) = ...
                    summarize(dat.(['z_prev' lname  'resp__' fullmodelname]) - ...
                    dat.(['z_prev' lname  'stim__' fullmodelname]));
                    end
                    try
                alldata.z_error(d,l) = ...
                    summarize(dat.(['z_prev' lname  'resp__' bestmodelname]) - ...
                    dat.(['z_prev' lname  'stim__' bestmodelname]));

                alldata.([vars{v} '_pval'])(d,l) = posteriorpval(traces.(['z_prev' lname  'resp']) - ...
                    traces.(['z_prev' lname  'stim']), 0);
                    end

                case 'v_correct'
                    try
                alldata.([vars{v} '_fullmodel'])(d,l) = ...
                    summarize(dat.(['v_prev' lname  'resp__' fullmodelname]) + ...
                    dat.(['v_prev' lname  'stim__' fullmodelname]));
                    end
                    try
                alldata.v_correct(d,l) = ...
                    summarize(dat.(['v_prev' lname  'resp__' bestmodelname]) + ...
                    dat.(['v_prev' lname  'stim__' bestmodelname]));

                alldata.([vars{v} '_pval'])(d,l) = posteriorpval(traces.(['v_prev' lname  'resp']) + ...
                    traces.(['v_prev' lname  'stim']), 0);
                    end

                case 'v_error'
                    try
                alldata.([vars{v} '_fullmodel'])(d,l) = ...
                    summarize(dat.(['v_prev' lname  'resp__' fullmodelname]) - ...
                    dat.(['v_prev' lname  'stim__' fullmodelname]));
                    end
                    try
                alldata.v_error(d,l) = ...
                    summarize(dat.(['v_prev' lname  'resp__' bestmodelname]) - ...
                    dat.(['v_prev' lname  'stim__' bestmodelname]));

                alldata.([vars{v} '_pval'])(d,l) = posteriorpval(traces.(['v_prev' lname  'resp']) - ...
                    traces.(['v_prev' lname  'stim']), 0);
                    end

            case 'v_prevresp'
                    try
                alldata.([vars{v} '_fullmodel'])(d,l) = ...
                    summarize(dat.(['v_prev' lname  'resp__' fullmodelname]));
                    end
                    try
                alldata.([vars{v}])(d,l) = ...
                    summarize(dat.(['v_prev' lname  'resp__' bestmodelname]));
                alldata.([vars{v} '_pval'])(d,l) = posteriorpval(traces.(['v_prev' lname  'resp']), 0);
                    end   

            case 'z_prevresp'
                    try
                alldata.([vars{v} '_fullmodel'])(d,l) = ...
                    summarize(dat.(['z_prev' lname  'resp__' fullmodelname]));
                    end
                    try
                alldata.([vars{v}])(d,l) = ...
                    summarize(dat.(['z_prev' lname  'resp__' bestmodelname]));
                alldata.([vars{v} '_pval'])(d,l) = posteriorpval(traces.(['z_prev' lname  'resp']), 0);
                    end   

            case 'v_prevstim'
                    try
                alldata.([vars{v} '_fullmodel'])(d,l) = ...
                    summarize(dat.(['v_prev' lname  'stim__' fullmodelname]));
                    end
                    try
                alldata.([vars{v}])(d,l)= ...
                    summarize(dat.(['v_prev' lname  'stim__' bestmodelname]));
                alldata.([vars{v} '_pval'])(d,l) = posteriorpval(traces.(['v_prev' lname  'stim']), 0);
                    end   

            case 'z_prevstim'
                    try
                alldata.([vars{v} '_fullmodel'])(d,l) = ...
                    summarize(dat.(['z_prev' lname  'stim__' fullmodelname]));
                    end
                    try
                alldata.([vars{v}])(d,l) = ...
                    summarize(dat.(['z_prev' lname  'stim__' bestmodelname]));
                alldata.([vars{v} '_pval'])(d,l) = posteriorpval(traces.(['z_prev' lname  'stim']), 0);
                    end   

            end % switch case

        end
    
    end
end

% ========================================================== %
% 3. PLOT THE VARIABLES THAT ARE PRESENT FOR THIS BEST MODEL
% ========================================================== %

% plot the thin lines only for weight that are not already in the bestmodel
for v = 1:length(vars),
    alldata.([vars{v} '_fullmodel'])(~isnan(alldata.([vars{v}]))) = ...
    alldata.([vars{v}])(~isnan(alldata.([vars{v}])));
end

colors = cbrewer('qual', 'Set2', length(datasets));

% CREATE FIGURE
for pltidx = 1:length(vars),
    
    close all;
    sp1 = subplot(4,4,1); hold on;
    plot([1 numlags], [0 0], 'k', 'linewidth', 0.5);
    
    for d = 1:length(datasets),
        % full model beneath, thin line
        plot(1:numlags, alldata.([vars{pltidx} '_fullmodel'])(d, :), 'color', colors(d, :), 'linewidth', 0.5);
        plot(1:numlags, alldata.(vars{pltidx})(d, :), 'color', colors(d, :), 'linewidth', 1);
    
        h = (alldata.([vars{pltidx} '_pval'])(d,:) < 0.05);
        if any(h>0),
           % plot(find(h==1), alldata.(vars{pltidx})(d, (h==1)), '.', 'markeredgecolor', colors(d, :), ...
           %    'markerfacecolor', colors(d,:), 'markersize', 7);
        end
    end

    % average across datasets
    plot(1:numlags, nanmean(alldata.([vars{pltidx} '_fullmodel'])), 'k', 'linewidth', 1);
    [h, pval] = ttest(alldata.([vars{pltidx} '_fullmodel']));
    [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pval);

    if any(adj_p < 0.05),
        plot(find(adj_p < 0.05), nanmean(alldata.([vars{pltidx} '_fullmodel'])(:, (adj_p < 0.05))), ...
            'k.', 'markersize', 10);
    end
    xlabel('Lags (# trials)');
    ylabel(regexprep(regexprep(regexprep(regexprep(vars{pltidx}, '_', ' ~ previous '), ...
        'v ', 'v_{bias} '), 'prevresp', 'response'), 'prevstim', 'stimulus'));
    set(gca, 'xtick', 1:numlags, 'xticklabel', lagnames, 'xcolor', 'k', 'ycolor', 'k');
    axis tight; offsetAxes;
    
    tightfig;
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/regressionkernels_correcterror_%d.pdf', pltidx));
    fprintf('~/Data/serialHDDM/regressionkernels_correcterror_%d.pdf \n', pltidx)
end


% ========================================================== %
% 4. PLOT CORRELATION KERNELS
% ========================================================== %

repMets = {'repetition', 'repetition_corrected', 'repetition-trivial', 'logistic', 'logistic-orth'};
for r = 1:length(repMets),

    disp(repMets{r})
    mat_z.r     = nan(length(datasets), numlags);
    mat_z.r_fullmodel     = nan(length(datasets), numlags);
    mat_z.pval  = nan(length(datasets), numlags);
    mat_dc.r    = nan(length(datasets), numlags);
    mat_dc.r_fullmodel    = nan(length(datasets), numlags);
    mat_dc.pval = nan(length(datasets), numlags);

    for d = 1:length(datasets),

        % from which model to take the correlations?
        bestmodelname = bestmodelnames{d};
        dat = readtable(sprintf('%s/summary/%s/allindividualresults.csv', mypath, datasets{d}));
        dat = dat(dat.session == 0, :);

        % covariates = zeros(size(dat.(['repetition' num2str(l)])));
        for l = 1:numlags,

            if l == 1,
                lname = [];
            elseif l == 7,
                lname = '7_10';
            elseif l == 8,
                lname = '11_15';
            else
                lname = num2str(l);
            end

            % what metric for repetition do we want to use?
            switch repMets{r}
            case 'repetition'
                repeat = dat.(['repetition' num2str(l)]);
            case 'repetition_corrected'
                repeat = dat.(['repetition_corrected' num2str(l)]);
            case 'repetition-trivial'
                repeat = dat.(['repetition' num2str(l)]) - arrayfun(@bernouilli_process, ...
                l*ones(size(dat.repetition)), dat.repetition);
                if l == 1,
                    repeat = dat.repetition;
                end
            case 'logistic'
                repeat = dat.(['logistic' num2str(l)]);
            case 'logistic-orth'
                repeat = dat.(['logistic_orth' num2str(l)]);
            end
            
            try
            % compute correlation coefficients
            [mat_z.r(d, l), mat_z.ci(d,l,:), mat_z.pval(d,l)] = ...
                spearmans(dat.(['z_prev' num2str(lname) 'resp__' bestmodelname]), ...
                repeat);

            [mat_dc.r(d, l), mat_dc.ci(d,l,:), mat_dc.pval(d,l)] = ...
                spearmans(dat.(['v_prev' num2str(lname) 'resp__' bestmodelname]), ...
                repeat);
            end
        end
    end

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
    set(gca, 'xtick', 1:numlags, 'xticklabel', lagnames, 'xcolor', 'k', 'ycolor', 'k');
    xlabel('Lags (# trials)');
    axis tight;
    ylim([-0.5 1]);
    offsetAxes; tightfig;
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/correlationkernels_z_%s.pdf', repMets{r}));

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
    set(gca, 'xtick', 1:numlags, 'xticklabel', lagnames, 'xcolor', 'k', 'ycolor', 'k');
    xlabel('Lags (# trials)');
    axis tight;
    ylim([-0.5 1]);
    offsetAxes; tightfig;
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/correlationkernels_dc_%s.pdf', repMets{r}));
    end

end

function y = summarize(x)

global individualrep

% flip weights around for alternators
y = nanmean(individualrep .* x);
% y = nanmean(x);

end