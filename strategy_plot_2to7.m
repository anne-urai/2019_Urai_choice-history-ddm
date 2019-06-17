function strategy_plot_2to7

% Code to fit the history-dependent drift diffusion models as described in
% Urai AE, de Gee JW, Tsetsos K, Donner TH (2019) Choice history biases subsequent evidence accumulation. eLife, in press.
%
% MIT License
% Copyright (c) Anne Urai, 2019
% anne.urai@gmail.com

global mypath datasets datasetnames
addpath(genpath('~/code/Tools'));
warning off; close all;

markers = {'d', 's', '^', 'v',  '>', '<'};
colors = cbrewer('qual', 'Set2', length(datasets));

for d = 1:length(datasets),
	 filename = dir(sprintf('%s/%s/*.csv', mypath, datasets{d}));
	 data  = readtable(sprintf('%s/%s/%s', mypath, datasets{d}, filename.name));
	 %assert(1==0)
	 if ~ismember(data.Properties.VariableNames, 'correct'),
	 	data.correct = (sign(data.stimulus) == data.response);
	 end

	avg_prevcorrect = nan(numel(unique(data.subj_idx)), 7);
	avg_preverror   = nan(numel(unique(data.subj_idx)), 7);

	for lag = 2:7,
		repeat 	= double(data.response == circshift(data.response, lag));
		% exclude trials that are not continuous
		trialnum = (data.trial == circshift(data.trial, lag)+lag);
		repeat(~trialnum) = NaN;

		% split into correct and error
		repeat_correct = repeat;
		repeat_correct(circshift(data.correct, lag) == 0) = NaN;
		repeat_error   = repeat;
		repeat_error(circshift(data.correct, lag) == 1) = NaN;

		avg_prevcorrect(:, lag) = splitapply(@nanmean, repeat_correct, findgroups(data.subj_idx));
		avg_preverror(:, lag) = splitapply(@nanmean, repeat_error, findgroups(data.subj_idx));
	end

	close all; subplot(441); hold on;
    plot([0.5 0.5], [0.3 0.7], 'color', 'k', 'linewidth', 0.5);
    plot([0.3 0.7], [0.5 0.5], 'color', 'k', 'linewidth', 0.5);
    scatter(nanmean(avg_prevcorrect, 2), nanmean(avg_preverror, 2), 5, colors(d, :), markers{d});
    xlabel({'P(repeat) after correct', 'lag 2-7'}); 
    ylabel({'P(repeat) after error', 'lag 2-7'});
    set(gca, 'xtick', 0.3:0.2:0.7, 'ytick', 0.3:0.2:0.7, 'xcolor', 'k', 'ycolor', 'k');
    xlim([0.3 0.7]); ylim([0.3 0.7]);
    offsetAxes; axis square;
    tightfig;
    if d < length(datasets), xlabel(''); end
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/strategyPlot_2to7_%d.pdf', d));

end

end