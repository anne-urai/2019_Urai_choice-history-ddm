function history_kernels

% Code to fit the history-dependent drift diffusion models as described in
% Urai AE, de Gee JW, Tsetsos K, Donner TH (2019) Choice history biases subsequent evidence accumulation. eLife, in press.
%
% MIT License
% Copyright (c) Anne Urai, 2019
% anne.urai@gmail.com

global mypath datasets datasetnames
addpath(genpath('~/code/Tools'));
warning off; close all;

%markers = {'d', 's', '^', 'v',  '>', '<'};
colors = cbrewer('qual', 'Set2', length(datasets));

for d = 1:length(datasets),
	 filename = dir(sprintf('%s/%s/*.csv', mypath, datasets{d}));
	 data  = readtable(sprintf('%s/%s/%s', mypath, datasets{d}, filename.name));

	avg = nan(numel(unique(data.subj_idx)), 7);
	for lag = 1:7,
		repeat 	= double(data.response == circshift(data.response, lag));
		% exclude trials that are not continuous
		trialnum = (data.trial == circshift(data.trial, lag)+lag);
		repeat(~trialnum) = NaN;

		avg(:, lag) = splitapply(@nanmean, repeat, findgroups(data.subj_idx));
	end

	close all; subplot(441); hold on;
	plot([0.5 7], [0.5 0.5], 'k');
	plot(1:7, avg, '-', 'color', [0.7 0.7 0.7]);
	plot(1:7, mean(avg), '-', 'color', colors(d, :), 'linewidth', 2);
	set(gca, 'xtick', 1:7);

    ylabel('P(repeat)')   
    if d == length(datasets),
        xlabel('Lag (# trials)');
    else
        set(gca, 'xticklabel', []);
    end
    %title(datasetnames{d});
    axis tight; offsetAxes;
    set(gca, 'xcolor', 'k', 'ycolor', 'k');
    tightfig;
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/historyKernel_d%d.pdf',d));

end

end