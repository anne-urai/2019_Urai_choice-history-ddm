function cumulative_prepeat()

addpath(genpath('~/code/Tools'));
warning off; % close all;
global datasets datasetnames mypath colors

alldat.alt = nan(6, 11);
alldat.rep = nan(6, 11);

alldat_alt.alt = nan(6, 11);
alldat_alt.rep = nan(6, 11);

for d = 1:length(datasets),

    dat = readtable(sprintf('%s/summary/%s/allindividualresults.csv', mypath, datasets{d}));
    dat = dat(dat.session == 0, :);

    % FOR ALTERNATORS, FLIP AROUND??

    cumrep = [dat.cum_rep0 dat.cum_rep1 dat.cum_rep2 dat.cum_rep3 dat.cum_rep4 dat.cum_rep5 ...
    dat.cum_rep6 dat.cum_rep7 dat.cum_rep8 dat.cum_rep9 dat.cum_rep10];

	cumalt = [dat.cum_alt0 dat.cum_alt1 dat.cum_alt2 dat.cum_alt3 dat.cum_alt4 dat.cum_alt5 dat.cum_alt6 ...
	dat.cum_alt7 dat.cum_alt8 dat.cum_alt9 dat.cum_alt10];

    % flip around for alternators!
    %cumrep(dat.repetition < 0.5, :) =  -1 * (cumrep(dat.repetition < 0.5, :) - 0.5) + 0.5;

    % split for repeaters and alternators
    alldat.alt(d, :) = nanmean(cumrep(dat.repetition < 0.5, :));
    alldat.rep(d, :) = nanmean(cumrep(dat.repetition > 0.5, :));
    alldat.all(d, :) = nanmean(cumrep);

    alldat_alt.alt(d, :) = nanmean(cumalt(dat.repetition < 0.5, :));
    alldat_alt.rep(d, :) = nanmean(cumalt(dat.repetition > 0.5, :));
    alldat_alt.all(d, :) = nanmean(cumalt);

end

colors = cbrewer('qual', 'Set2', length(datasets));
groups = {'all', 'rep', 'alt'};

close all;
for g = 1:length(groups),
subplot(4,4,g); hold on;

    plot([0, 10], [0.5 0.5], 'k', 'linewidth', 0.5);

    for d = 1:length(datasets),
        % full model beneath, thin line
        plot(0:10, alldat.(groups{g})(d, :), 'color', colors(d, :), 'linewidth', 0.5);
    end

ylim([0.4 0.7]);
errorbar(0:10, nanmean(alldat.(groups{g})), nanstd(alldat.(groups{g})) ./ sqrt(6), ...
	 'ok-', 'linewidth', 1, 'capsize', 0, 'markersize', 3, 'markeredgecolor', 'k', 'markerfacecolor', 'k');

% add stats
[h, pval] = ttest(alldat.(groups{g}) - 0.5);
[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pval);
if any(adj_p < 0.05),
    plot(find(adj_p < 0.05)-1, 0.45*ones(length(find(adj_p < 0.05)), 1), ...
        'k.', 'markersize', 10);
end

switch groups{g}
    case 'all'
                title('All observers', 'fontweight', 'normal', 'fontangle', 'italic');
	case 'rep'
		title('Repeaters', 'fontweight', 'normal', 'fontangle', 'italic');
	case 'alt'
	title('Alternators', 'fontweight', 'normal', 'fontangle', 'italic');
end
	
set(gca, 'xtick', 0:6);
xlim([0 6]);
offsetAxes;
xlabel('Preceding # repetitions');
if g == 1, ylabel('P(repeat)'); end
end

tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/cumulative_prepeat.pdf'));


%% PRINT AGAIN, NOW FOR ALTERNATING SEQUENCES

close all;
for g = 1:length(groups),
subplot(4,4,g); hold on;

    plot([0, 10], [0.5 0.5], 'k', 'linewidth', 0.5);

    for d = 1:length(datasets),
        % full model beneath, thin line
        plot(0:10, alldat_alt.(groups{g})(d, :), 'color', colors(d, :), 'linewidth', 0.5);
    end

errorbar(0:10, nanmean(alldat_alt.(groups{g})), nanstd(alldat_alt.(groups{g})) ./ sqrt(6), ...
	 'ok-', 'linewidth', 1, 'capsize', 0, 'markersize', 3, 'markeredgecolor', 'k', 'markerfacecolor', 'k');

ylim([0.35 0.6]);
% assert(1==0)
% add stats
[h, pval] = ttest(alldat_alt.(groups{g}) - 0.5);
[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pval);
if any(adj_p < 0.05),
    plot(find(adj_p < 0.05)-1, 0.57*ones(length(find(adj_p < 0.05)), 1), ...
        'k.', 'markersize', 10);
end

switch groups{g}
case 'all'
                    title('All observers', 'fontweight', 'normal', 'fontangle', 'italic');
	case 'rep'
		title('Repeaters', 'fontweight', 'normal', 'fontangle', 'italic');
	case 'alt'
	title('Alternators', 'fontweight', 'normal', 'fontangle', 'italic');
end
	
xlim([0 6]);
set(gca, 'xtick', 0:6);
offsetAxes;
xlabel('Preceding # alternations');
if g == 1, ylabel('P(repeat)'); end

end
tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/cumulative_prepeat_alternation.pdf'));

    %%%%%%%%%%%%%%%%%%%%%%%%

close all;
for g = 1:length(groups),
subplot(4,4,g); hold on;

    plot([0, 10], [0.5 0.5], 'k', 'linewidth', 0.5);

    errorbar(0:10, nanmean(alldat.(groups{g})), nanstd(alldat.(groups{g})) ./ sqrt(6), ...
     'ob-', 'linewidth', 1, 'capsize', 0, 'markersize', 3, 'markeredgecolor', 'b', 'markerfacecolor', 'b');

[h, pval] = ttest(alldat.(groups{g}) - 0.5);
[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pval);
if any(adj_p < 0.05),
    plot(find(adj_p < 0.05)-1, 0.65*ones(length(find(adj_p < 0.05)), 1), ...
        '.b',  'markersize', 10);
end

errorbar(0:10, nanmean(alldat_alt.(groups{g})), nanstd(alldat_alt.(groups{g})) ./ sqrt(6), ...
         'or-', 'linewidth', 1, 'capsize', 0, 'markersize', 3, 'markeredgecolor', 'r', 'markerfacecolor', 'r');

[h, pval] = ttest(alldat_alt.(groups{g}) - 0.5);
[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pval);
if any(adj_p < 0.05),
    plot(find(adj_p < 0.05)-1, 0.4*ones(length(find(adj_p < 0.05)), 1), ...
        'r.', 'markersize', 10);
end

% lag 1 in black
errorbar(0, nanmean(alldat_alt.(groups{g})(:, 1)), nanstd(alldat_alt.(groups{g})(:, 1), [], 1) ./ sqrt(6), ...
         'ok-', 'linewidth', 1, 'capsize', 0, 'markersize', 3, 'markeredgecolor', 'k', 'markerfacecolor', 'k');

switch groups{g}
case 'all'
     title('All observers', 'fontweight', 'normal', 'fontangle', 'italic');
    case 'rep'
        title('Repeaters', 'fontweight', 'normal', 'fontangle', 'italic');
    case 'alt'
    title('Alternators', 'fontweight', 'normal', 'fontangle', 'italic');
end
    
xlim([0 6]); ylim([0.4 0.7]);
set(gca, 'xtick', 0:6);
offsetAxes;
xlabel('Preceding sequence length');
if g == 1, ylabel('P(repeat)'); end

end
tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/cumulative_prepeat_combined.pdf'));





end