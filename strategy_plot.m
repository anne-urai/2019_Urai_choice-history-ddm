function strategyPlot
% Code to fit the history-dependent drift diffusion models described in
% Urai AE, Gee JW de, Donner TH (2018) Choice history biases subsequent evidence accumulation. bioRxiv:251595
%
% MIT License
% Copyright (c) Anne Urai, 2018
% anne.urai@gmail.com

global mypath datasets datasetnames
addpath(genpath('~/code/Tools'));
warning off; close all;

cmap = viridis(256);
colormap(cmap);


for d = length(datasets):-1:1,
    close all;
    subplot(441); hold on;
    markers = {'d', 's', '^', 'v',  '>', '<'};
    colors = cbrewer('qual', 'Set2', length(datasets));
    
    plot([0.5 0.5], [0.3 0.7], 'color', 'k', 'linewidth', 0.5);
    plot([0.3 0.7], [0.5 0.5], 'color', 'k', 'linewidth', 0.5);
    
    dat = readtable(sprintf('%s/summary/%s/allindividualresults.csv', mypath, datasets{d}));
    dat = dat(dat.session == 0, :);
    s{d} = scatter(dat.repetition_prevcorrect, dat.repetition_preverror, 5, colors(d, :), markers{d});
    legtxt{d} = cat(2, datasetnames{d}{1}, ' ', datasetnames{d}{2});
    alldat{d} = dat(:, {'repetition_prevcorrect', 'repetition_preverror'});
    
    xlabel('P(repeat) after correct');
    ylabel('P(repeat) after error');
    set(gca, 'xtick', 0.3:0.2:0.7, 'ytick', 0.3:0.2:0.7, 'xcolor', 'k', 'ycolor', 'k');
    xlim([0.3 0.7]); ylim([0.3 0.7]);
    offsetAxes; axis square;
    
    tightfig;
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/strategyPlot_%d.pdf', d));
end

% display how many observers are in each quadrant!
alldat = cat(1, alldat{:});

stay = sum(alldat.repetition_prevcorrect > 0.5 & alldat.repetition_preverror > 0.5)
alternate = sum((alldat.repetition_prevcorrect < 0.5 & alldat.repetition_preverror < 0.5))
winstay_loseswitch = sum(alldat.repetition_prevcorrect > 0.5 & alldat.repetition_preverror < 0.5)
winswitch_losestay = sum(alldat.repetition_prevcorrect < 0.5 & alldat.repetition_preverror > 0.5)

end
