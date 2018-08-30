% SIMULATE BETWEEN-SUBJECT CORRELATIONS BETWEEN TWO VARIABLES AS A FUNCTION
% OF THE NUMBER OF TRIALS AT THE INDIVIDUAL LEVEL %
%
% RELATED TO FIGURE 5: IS THE CORRELATION FOR ERROR TRIALS LOWER BECAUSE OF
% THE TRIAL COUNT?

% ============================================================== %
% GENERATE SUBJECTS WITH EACH AN INDIVIDUAL BEHAVIORAL TENDENCY
% ============================================================== %

% annotate with Anke's dataset
data = readtable('~/Data/HDDM/Anke_MEG_transition/Anke_MEG_transition.csv');
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
    
    for n = 1:length(ntrials),
        
        % simulate two variables that are independently generated based on the same
        % generative probability
        var1 = binornd(ntrials(n), true_rep_prob) ./ ntrials(n);
        
        % this is key: how much noise is added in the estimation of
        % history-shift in drift bias?
        var2 = binornd(ntrials(n), true_rep_prob) ./ ntrials(n);
        
        between_sj_corr_pearson(n)  = corr(var1, var2);
        between_sj_corr_spearman(n) = corr(var1, var2, 'type', 'spearman');
        
    end
    
    % ============================================================== %
    % PLOT
    % ============================================================== %
    
    subplot(2,2,sj);
    p1 = plot(ntrials, between_sj_corr_pearson);
    hold on;
    grid on;
    p2 = plot(ntrials, between_sj_corr_spearman);
    xlabel('# trials per individual');
    ylabel({'Correlation coefficient', sprintf('across %d subjects', sjgroups(sj))});
    box off; ylim([0.8 1]); set(gca, 'ytick', [0.8:0.025:1]);
    legend([p1 p2], {'Pearson', 'Spearman'}, 'box', 'off', 'location', 'southeast', 'AutoUpdate','off');
    
    if sj == 1,
        % add
        vline(mean(trialcount(correct == 0)), 'color', 'k');
        text(mean(trialcount(correct == 0)), 1.01, 'error trial count', ...
            'horizontalalignment', 'center', 'backgroundcolor', 'w', 'fontsize', 6);
        vline(mean(trialcount(correct == 1)), 'color', 'k');
        text(mean(trialcount(correct == 1)), 1.01, 'correct trial count', ...
            'horizontalalignment', 'center', 'backgroundcolor', 'w', 'fontsize', 6);
    end
    
end
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/trialcount_correlation_simulation.pdf'));

