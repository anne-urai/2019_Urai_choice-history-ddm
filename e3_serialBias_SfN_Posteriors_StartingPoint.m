function e3_serialBias_SfN_Posteriors_StartingPoint

% Code to fit the history-dependent drift diffusion models described in
% Urai AE, Gee JW de, Donner TH (2018) Choice history biases subsequent evidence accumulation. bioRxiv:251595
%
% MIT License
% Copyright (c) Anne Urai, 2018
% anne.urai@gmail.com

addpath(genpath('~/code/Tools'));
warning off; close all;
global mypath datasets datasetnames colors

colors2 = cbrewer('qual', 'Set2', length(datasets));

% STARTING POINT POSTERIORS 
close; subplot(3,3,1); hold on;
plot([1 6], [0 0], 'k-', 'linewidth', 0.5);

for d = 1:length(datasets),
    dat = readtable(sprintf('%s/%s/stimcoding_dc_z_prevresp/group_traces.csv', mypath, datasets{d}));
    difference = dat.z_trans_1_ - dat.z_trans__1_;
    h = violinPlot(difference, 'color', colors2(d, :), 'showMM', 6, 'xValues', d);
    legtext{d} = cat(2, datasetnames{d}{1}, ' ', datasetnames{d}{2});
end

set(gca, 'xtick', 1:length(datasets), 'xticklabel', legtext, 'xticklabelrotation', -30, 'xcolor', 'k');
ylabel({'History shift in z'}, 'color', colors(1, :));
offsetAxes;
tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/bias_posteriors_startingpoint.pdf'));


% DRIFT BIAS POSTERIORS
close; subplot(3,3,1); hold on;
plot([1 6], [0 0], 'k-', 'linewidth', 0.5);

for d = 1:length(datasets),
    dat = readtable(sprintf('%s/%s/stimcoding_dc_z_prevresp/group_traces.csv', mypath, datasets{d}));
    difference = dat.dc_1_ - dat.dc__1_;
    h = violinPlot(difference, 'color', colors2(d, :), 'showMM', 6, 'xValues', d);
    legtext{d} = cat(2, datasetnames{d}{1}, ' ', datasetnames{d}{2});
end

set(gca, 'xtick', 1:length(datasets), 'xticklabel', legtext, 'xticklabelrotation', -30, 'xcolor', 'k');
ylabel({'History shift in v_{bias}'}, 'color', colors(2, :));
offsetAxes;
tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/bias_posteriors_drift bias.pdf'));

end
