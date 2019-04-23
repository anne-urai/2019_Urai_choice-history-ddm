function kernels_lags_bestmodel

global mypath datasets 
addpath(genpath('~/code/Tools'));
warning off; close all;

numlags = 6;
vars = {'z_correct', 'z_error', 'v_correct', 'v_error', ...
    'z_prevresp', 'z_prevstim', 'v_prevresp', 'v_prevstim'};
lagnames = {'1', '2', '3', '4', '5', '6'};
fullmodelname = 'regressdczlag6'; % extend thin lines for weights from biggest model
global individualrep

for m = 1:length(vars),
    alldata.(vars{m})       = nan(length(datasets), numlags);
    alldata.([vars{m} '_fullmodel'])       = nan(length(datasets), numlags);
    alldata.([vars{m} '_pval'])   = nan(length(datasets), numlags);
end

for d = 1:length(datasets),

    % ========================================================== %
    % GET WEIGHTS FOR THE FULL MODEL 
    % ========================================================== %
    
    dat = readtable(sprintf('%s/summary/%s/allindividualresults.csv', mypath, datasets{d}));
    dat = dat(dat.session == 0, :);

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

assert(1==0)

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
    set(gca, 'xtick', 1:numlags, 'xticklabel', lagnames, 'xcolor', 'k', 'ycolor', 'k');
    axis tight; offsetAxes;
    
    tightfig;
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/regressionkernels_correcterror_%d.pdf', pltidx));
    % fprintf('~/Data/serialHDDM/regressionkernels_correcterror_%d.pdf \n', pltidx)
end


end

function y = summarize(x)

global individualrep

% flip weights around for alternators
y = nanmean(individualrep .* x);
% y = nanmean(x);

end