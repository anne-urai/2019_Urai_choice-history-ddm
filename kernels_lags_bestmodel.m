function kernels_lags_bestmodel

global mypath datasets
addpath(genpath('~/code/Tools'));
warning off; close all;


numlags = 6;
lagnames = {'1', '2', '3', '4', '5', '6', '7-10', '11-15'};
vars = {'z_correct', 'z_error', 'v_correct', 'v_error', ...
    'z_prevresp', 'z_prevstim', 'v_prevresp', 'v_prevstim'};
for m = 1:length(vars),
    alldata.(vars{m})                   = nan(length(datasets), numlags);
    alldata.([vars{m} '_fullmodel'])    = nan(length(datasets), numlags);
    alldata.([vars{m} '_pval'])         = nan(length(datasets), numlags);
    alldata.([vars{m} '_ylabel'])       = cell(length(datasets));
end
fullmodelname = 'regressdczlag6'; % extend thin lines for weights from biggest model
global individualrep

global flipAlternators;
flipAlternators = 1;

for d = 1:length(datasets),
    
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
    
    % ============================= %
    % 1. DETERMINE THE BEST MODEL
    % ============================= %
    
    mdldic = nan(1, length(mdls));
    for m = 1:length(mdls),
        try
            modelcomp = readtable(sprintf('%s/%s/%s/model_comparison.csv', ...
                mypath, datasets{d}, mdls{m}), 'readrownames', true);
            mdldic(m) = modelcomp.aic;
        catch
            fprintf('%s/%s/%s/model_comparison.csv  NOT FOUND\n', ...
                mypath, datasets{d}, mdls{m})
        end
    end
    
    % everything relative to the full model
    mdldic = bsxfun(@minus, mdldic, mdldic(1));
    mdldic = mdldic(2:end);
    mdls = mdls(2:end);
    [~, bestMdl] = min(mdldic);
    
    % everything relative to the full model
    bestmodelname = regexprep(regexprep(mdls{bestMdl}, '_', ''), '-', 'to');
    bestmodelnames{d} = bestmodelname;
    disp(bestmodelname);
    
    % ========================================================== %
    % 2. FOR THIS MODEL, RECODE INTO CORRECT AND ERROR
    % ========================================================== %
    
    dat = readtable(sprintf('%s/summary/%s/allindividualresults.csv', mypath, datasets{d}));
    dat = dat(dat.session == 0, :);
    try
        traces = readtable(sprintf('%s/%s/%s/group_traces.csv', mypath, datasets{d}, mdls{bestMdl+1}));
    end
    % flip around weights for alternators
    individualrep = sign(dat.repetition - 0.5);
    
    for l = 1:numlags,
        if l == 1,
            lname = '';
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
                    
                    if flipAlternators,
                        alldata.([vars{v} '_ylabel']){d} = {'z ~ history bias', 'previous correct'};
                    else
                        alldata.([vars{v} '_ylabel']){d} = {'z ~ repetition', 'previous correct'};
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
                    if flipAlternators,
                        alldata.([vars{v} '_ylabel']){d} = {'z ~ history bias', 'previous error'};
                    else
                        alldata.([vars{v} '_ylabel']){d} = {'z ~ repetition', 'previous error'};

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
                    
                    if flipAlternators,
                        alldata.([vars{v} '_ylabel']){d} = {'v_{bias} ~ history bias', 'previous correct'};
                    else
                        alldata.([vars{v} '_ylabel']){d} = {'v_{bias} ~ repetition', 'previous correct'};
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
                    
                    if flipAlternators,
                        alldata.([vars{v} '_ylabel']){d} = {'v_{bias} ~ history bias', 'previous error'};
                    else 
                        alldata.([vars{v} '_ylabel']){d} = {'v_{bias} ~ repetition', 'previous error'};
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
                    
                    if flipAlternators,
                        alldata.([vars{v} '_ylabel']){d} = 'v_{bias} ~ history bias';
                    else
                        alldata.([vars{v} '_ylabel']){d} = 'v_{bias} ~ repetition';
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
                    
                    if flipAlternators,
                        alldata.([vars{v} '_ylabel']){d} = 'z ~ history bias';
                    else
                        alldata.([vars{v} '_ylabel']){d} = 'z ~ repetition';
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
                    
                    if flipAlternators,
                        alldata.([vars{v} '_ylabel']){d} = 'v_{bias} ~ history bias in previous stim';
                    else
                        alldata.([vars{v} '_ylabel']){d} = 'v_{bias} ~ previous stim';
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
                    
                    
                    if flipAlternators,
                        alldata.([vars{v} '_ylabel']){d} = 'z ~ history bias in previous stim';
                    else
                        alldata.([vars{v} '_ylabel']){d} = 'z ~ previous stim';
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
        plot(1:numlags, alldata.([vars{pltidx} '_fullmodel'])(d, :), 'color', colors(d, :), 'linewidth', 0.2);
        plot(1:numlags, alldata.(vars{pltidx})(d, :), 'color', colors(d, :), 'linewidth', 1);
        
        % h = (alldata.([vars{pltidx} '_pval'])(d,:) < 0.05);
        % if any(h>0),
        %    % plot(find(h==1), alldata.(vars{pltidx})(d, (h==1)), '.', 'markeredgecolor', colors(d, :), ...
        %    %    'markerfacecolor', colors(d,:), 'markersize', 7);
        % end
    end
    
    % average across datasets
    plot(1:numlags, nanmean(alldata.([vars{pltidx} '_fullmodel'])), 'k', 'linewidth', 1);
    % [h, adj_p] = ttest(alldata.([vars{pltidx}])); % stats on best fits
    % %[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pval);
    
    % if any(adj_p < 0.05),
    %     plot(find(adj_p < 0.05), nanmean(alldata.([vars{pltidx} '_fullmodel'])(:, (adj_p < 0.05))), ...
    %         'k.', 'markersize', 10);
    % end
    
    xlabel('Lags (# trials)');
    ylabel(regexprep(regexprep(regexprep(regexprep(vars{pltidx}, '_', ' ~ previous '), ...
        'v ', 'v_{bias} '), 'prevresp', 'response'), 'prevstim', 'stimulus'));
    ylabel(alldata.([vars{pltidx} '_ylabel']){d});
    
    set(gca, 'xtick', 1:numlags, 'xticklabel', lagnames, 'xcolor', 'k', 'ycolor', 'k');
    axis tight;
    ylim([-0.1 0.22]); offsetAxes;
    
    tightfig;
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/regressionkernels_correcterror_%d_flipAlt%d.pdf', pltidx, flipAlternators));
    % fprintf('~/Data/serialHDDM/regressionkernels_correcterror_%d.pdf \n', pltidx)
end


% ========================================================== %
% 4. REPEAT, BUT FIT EXPONENTIAL
% ========================================================== %

% CREATE FIGURE
for pltidx = 1:length(vars),
    
    
    %% FIT EXPONENTIAL
    f = fit(transpose(1:numlags), nanmean(alldata.([vars{pltidx} '_fullmodel']))','exp1');
    disp(f);
    
    % close all;
    % plot(f, transpose(1:numlags), nanmean(alldata.([vars{pltidx} '_fullmodel']))');
    % print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/regressionkernels_correcterror_exponential_%d_test.pdf', pltidx));
    
    close all;
    sp1 = subplot(4,4,1); hold on;
    plot([1 numlags], [0 0], 'k', 'linewidth', 0.5);
    
    for d = 1:length(datasets),
        % full model beneath, thin line
        plot(1:numlags, alldata.([vars{pltidx} '_fullmodel'])(d, :), 'color', colors(d, :), 'linewidth', 0.2);
        plot(1:numlags, alldata.(vars{pltidx})(d, :), 'color', colors(d, :), 'linewidth', 1);
    end
    
    
    errorbar(1:numlags, nanmean(alldata.([vars{pltidx} '_fullmodel'])), ...
        nanstd(alldata.([vars{pltidx} '_fullmodel'])) ./ sqrt(6), ...
        'ok', 'linewidth', 1, 'capsize', 0, 'markersize', 3, 'markeredgecolor', 'k', 'markerfacecolor', 'w');
    
    % INDICATE WHICH ARE SIGNIFICANT!
    [h, pval] = ttest(alldata.([vars{pltidx} '_fullmodel']));
    
    errorbar((pval < 0.05), nanmean(alldata.([vars{pltidx} '_fullmodel'])(:, pval < 0.05), ...
        nanstd(alldata.([vars{pltidx} '_fullmodel'])(:, pval < 0.05) ./ sqrt(6), ...
        'ok', 'linewidth', 1, 'capsize', 0, 'markersize', 3, 'markeredgecolor', 'w', 'markerfacecolor', 'k');
    
    % show the model on top
    params = coeffvalues(f);
    disp(coeffnames(f));
    y = feval(f, linspace(1, numlags, 100));
    plot(linspace(1, numlags, 100), y, 'color', 'k', 'linewidth', 1);
    text(3.3, -0.08, sprintf('$V(t) = %.2f e^{-t/%.2f}$', params(1), -1/params(2)), ...
        'interpreter', 'latex', 'fontsize', 5);
    %text(2, 0.08, sprintf('$\frac{%d}{%d}$', 1, 2),'Interpreter','latex')

    xlabel('Lags (# trials)');
    ylabel(regexprep(regexprep(regexprep(regexprep(vars{pltidx}, '_', ' ~ previous '), ...
        'v ', 'v_{bias} '), 'prevresp', 'response'), 'prevstim', 'stimulus'));
    ylabel(alldata.([vars{pltidx} '_ylabel']){d});
    
    
    set(gca, 'xtick', 1:numlags, 'xticklabel', lagnames, 'xcolor', 'k', 'ycolor', 'k');
    axis tight;
    ylim([-0.1 0.22]); offsetAxes;
    
    tightfig;
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/regressionkernels_correcterror_exponential_%d_flipAlt%d.pdf', pltidx, flipAlternators));
    
end
end

function y = summarize(x)

global individualrep
global flipAlternators

% flip weights around for alternators
if flipAlternators,
    y = nanmean(individualrep .* x);
else
    y = nanmean(x);
end

end
