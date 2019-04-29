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

    alldat_alt.alt(d, :) = nanmean(cumalt(dat.repetition < 0.5, :));
    alldat_alt.rep(d, :) = nanmean(cumalt(dat.repetition > 0.5, :));
end

colors = cbrewer('qual', 'Set2', length(datasets));

groups = {'rep', 'alt'};

close all;
for g = 1:length(groups),
subplot(4,4,g); hold on;

    plot([0, 10], [0.5 0.5], 'k', 'linewidth', 0.5);

    for d = 1:length(datasets),
        % full model beneath, thin line
        plot(0:10, alldat.(groups{g})(d, :), 'color', colors(d, :), 'linewidth', 0.5);
    end

errorbar(0:10, nanmean(alldat.(groups{g})), nanstd(alldat.(groups{g})) ./ sqrt(6), ...
	 'ok-', 'linewidth', 1, 'capsize', 0, 'markersize', 3, 'markeredgecolor', 'k', 'markerfacecolor', 'k');
axis tight;

switch groups{g}
	case 'rep'
		title('Repeaters', 'fontweight', 'normal', 'fontangle', 'italic');
	case 'alt'
	title('Alternators', 'fontweight', 'normal', 'fontangle', 'italic');
end
	
xlim([0 6]);
offsetAxes;
xlabel('Repeating sequence length');
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

errorbar(0:10, nanmean(alldat_alt.(groups{g})), nanstd(alldat.(groups{g})) ./ sqrt(6), ...
	 'ok-', 'linewidth', 1, 'capsize', 0, 'markersize', 3, 'markeredgecolor', 'k', 'markerfacecolor', 'k');
axis tight;

switch groups{g}
	case 'rep'
		title('Repeaters', 'fontweight', 'normal', 'fontangle', 'italic');
	case 'alt'
	title('Alternators', 'fontweight', 'normal', 'fontangle', 'italic');
end
	
xlim([0 6]);
offsetAxes;
xlabel('Alternating sequence length');
if g == 1, ylabel('P(repeat)'); end

end

tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/cumulative_prepeat_alternation.pdf'));


end