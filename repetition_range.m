function e3_serialBias_SfN_repetitionRange

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

for d = 1:length(datasets),
    
    close all;
    subplot(4,4,1); hold on;
    dat = readtable(sprintf('%s/summary/%s/allindividualresults.csv', mypath, datasets{d}));
    dat = dat(dat.session == 0, :);
    if d == 5
        dat(dat.subjnr == 11, :) = [];
    end
    
    rep = sort(dat.repetition);
    b = barh(rep, 'facecolor', [0.5 0.5 0.5], 'basevalue', 0.5, 'edgecolor', 'none');
    plot([0.5 0.5],[0.5 numel(rep)+0.5],  'k');
    b(1).BaseLine.LineStyle = 'none';
    
    % ylabel(sprintf('%s, n = %d', datasetnames{d}{1}, numel(rep)));
    % show on x-axis what the mean is
    box off; xlim([0.4 0.65]); ylim([1 numel(rep)+0.5]);
    plot(mean(rep), 0, 'k^', 'markerfacecolor', linspecer(1), 'markeredgecolor', linspecer(1), 'markersize', 3);
    set(gca, 'ytick', [1 numel(rep)], 'xtick', [0.4 0.5 0.6]);
    
    ylabel('# Observers')
    % % yyaxis right;
    % y = ylabel(datasetnames{d});s
    
    % y.Rotation = y.Rotation + 180;
    % y.Position(1) = y.Position(1) + 0.27;
    % y.Position(2) = y.Position(2) - 2;
    
    if d == length(datasets),
    xlabel('P(repeat)');
end
    offsetAxes;
    
    %axis square;
    set(gca, 'xcolor', 'k', 'ycolor', 'k');
    % text(0.7, 0, '.w');
    subplot(4,20,5); plot(0,0,'.w', 'color', 'w'); axis off;
    % offsetAxes;
    tightfig;
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/repetitionRange_d%d.pdf',d));
    % print(gcf, '-depsc', sprintf('~/Data/serialHDDM/repetitionRange_d%d.eps',d));
    
end