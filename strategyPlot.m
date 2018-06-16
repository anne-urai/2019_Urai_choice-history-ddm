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

close all;
subplot(3,3,1); hold on;
markers = {'d', 's', '^', 'v',  '>', '<'};
colors = cbrewer('qual', 'Set2', length(datasets));

plot([0.5 0.5], [0.3 0.7], 'color', 'k', 'linewidth', 0.5);
plot([0.3 0.7], [0.5 0.5], 'color', 'k', 'linewidth', 0.5);

for d = length(datasets):-1:1,
    dat = readtable(sprintf('%s/summary/%s/allindividualresults.csv', mypath, datasets{d}));
    dat = dat(dat.session == 0, :);
    s{d} = scatter(dat.repetition_prevcorrect, dat.repetition_preverror, 5, colors(d, :), markers{d});
    legtxt{d} = cat(2, datasetnames{d}{1}, ' ', datasetnames{d}{2});
end
xlabel('P(repeat) after correct');
ylabel('P(repeat) after error');
set(gca, 'xtick', 0.3:0.2:0.7, 'ytick', 0.3:0.2:0.7);
offsetAxes; axis square;
% 
% l = legend([s{:}], legtxt);
% l.Box = 'off';
% subplot(3,3,4); plot(1,1, '.w'); axis off;
% % hline(0.5); vline(0.5);
% tightfig;
% l.Position(2) = l.Position(2)  - .5;
% l.Position(1) = l.Position(1)  + 0.1;

tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/strategyPlot.pdf'));

end
