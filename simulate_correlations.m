% SIMULATE BETWEEN-SUBJECT CORRELATIONS BETWEEN TWO VARIABLES AS A FUNCTION
% OF THE NUMBER OF TRIALS AT THE INDIVIDUAL LEVEL %
%
% RELATED TO FIGURE 5: IS THE CORRELATION FOR ERROR TRIALS LOWER BECAUSE OF
% THE TRIAL COUNT?

% ============================================================== %
% GENERATE SUBJECTS WITH EACH AN INDIVIDUAL BEHAVIORAL TENDENCY
% ============================================================== %

% annotate with Anke's dataset
global mypath datasets
d = 2;

% load data
csvfile = dir(sprintf('%s/%s/*.csv', mypath, datasets{d}));
csvfile = csvfile(arrayfun(@(x) ~strcmp(x.name(1),'.'), csvfile)); % remove hidden files
data = readtable(sprintf('%s/%s/%s', mypath, datasets{d}, csvfile.name));


[gr, sj, correct] = findgroups(data.subj_idx, (sign(data.stimulus) == sign(data.response-0.5)));
trialcount = splitapply(@numel, data.response, gr);

% simulation
sjgroups = [32 100]; % try with different number of subjects
figure;
for sj =  1:length(sjgroups)
    true_rep_prob = normrnd(0.5, 0.1, sjgroups(sj), 1);
    
    % FOR EACH SUBJECT, GENERATE NTR TRIALS BASED ON THIS WEIGHTED COIN FLIP
    ntrials                     = 100:10:2000; % some steps
    between_sj_corr_pearson     = nan(1, length(ntrials));
    between_sj_corr_spearman    = nan(1, length(ntrials));
    between_sj_corr_slope       = nan(1, length(ntrials));
    
    for n = 1:length(ntrials),
        
        % simulate two variables that are independently generated based on the same
        % generative probability
        var1 = binornd(ntrials(n), true_rep_prob) ./ ntrials(n);
        
        % this is key: how much noise is added in the estimation of
        % history-shift in drift bias?
        var2 = binornd(ntrials(n), true_rep_prob) ./ ntrials(n);
        
        between_sj_corr_pearson(n)  = corr(var1, var2);
        between_sj_corr_spearman(n) = corr(var1, var2, 'type', 'spearman');
        p = polyfit(var1, var2, 1);
        between_sj_corr_slope(n) = p(1);
    end
    
    % ============================================================== %
    % PLOT
    % ============================================================== %
    
    subplot(2,2,sj);
    hold on;
    grid on;
    % plot(ntrials, between_sj_corr_slope);
    % plot(ntrials, between_sj_corr_pearson);
    plot(ntrials, between_sj_corr_spearman);
    xlabel('# trials per individual');
    ylabel({'Spearman''s\rho ', sprintf('across %d subjects', sjgroups(sj))});
    box off; ylim([0.8 1]); set(gca, 'ytick', [0.8:0.025:1]);
    % legend({'Regression slope', 'Pearson', 'Spearman'}, 'box', 'off', 'location', 'southeast', 'AutoUpdate','off');
    
    if sj == 1,
        % add
        vline(mean(trialcount(correct == 0)), 'color', 'k');
        text(mean(trialcount(correct == 0)), 0.85, 'error trial count', ...
            'horizontalalignment', 'center', 'fontsize', 7);
        vline(mean(trialcount(correct == 1)), 'color', 'k');
        text(mean(trialcount(correct == 1)), 0.9, 'correct trial count', ...
            'horizontalalignment', 'center','fontsize', 7);
    end
    
end
tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/trialcount_correlation_simulation.pdf'));

