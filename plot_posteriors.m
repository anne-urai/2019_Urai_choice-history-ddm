function plot_posteriors

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

% ========================================== %
% STARTING POINT POSTERIORS
% ========================================== %

close; subplot(3,3,1); hold on;
plot([1 6], [0 0], 'k-', 'linewidth', 0.5);

for d = 1:length(datasets),
    dat = readtable(sprintf('%s/%s/stimcoding_dc_z_prevresp/group_traces.csv', mypath, datasets{d}));
    difference = dat.z_trans_1_ - dat.z_trans__1_;
        violinPlot_distribution(d, difference, colors2(d, :));
    legtext{d} = cat(2, datasetnames{d}{1}, ' ', datasetnames{d}{2});
    pval = posteriorpval(dat.z_trans_1_, dat.z_trans__1_);

    if pval < 0.001,
            text(d+0.1, -0.5, sprintf('p < 0.001'), ...
        'fontangle', 'italic', 'fontsize', 4, 'rotation', -30);
        set(gca, 'xcolor', 'k', 'ycolor', 'k');
    else
    text(d+0.1, -0.5, sprintf('p = %.3f', pval), ...
        'fontangle', 'italic', 'fontsize', 4, 'rotation', -30);
        set(gca, 'xcolor', 'k', 'ycolor', 'k');
    end

end

set(gca, 'xtick', 1:length(datasets), 'xticklabel', legtext, 'xticklabelrotation', -30, 'xcolor', 'k');
ylabel({'History shift in z'}, 'color', colors(1, :));
axis tight; offsetAxes;
tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/bias_posteriors_startingpoint.pdf'));


% ========================================== %
% DRIFT BIAS POSTERIORS
% ========================================== %

close; subplot(3,3,1); hold on;
plot([1 6], [0 0], 'k-', 'linewidth', 0.5);

for d = 1:length(datasets),
    dat = readtable(sprintf('%s/%s/stimcoding_dc_z_prevresp/group_traces.csv', mypath, datasets{d}));
	disp(datasets{d});
    difference = dat.dc_1_ - dat.dc__1_;
    violinPlot_distribution(d, difference, colors2(d, :));
    legtext{d} = cat(2, datasetnames{d}{1}, ' ', datasetnames{d}{2});
    pval = posteriorpval(dat.dc_1_, dat.dc__1_);
    if pval < 0.001,
    text(d+0.1, -0.3, sprintf('p < 0.001'), ...
        'fontangle', 'italic', 'fontsize', 4, 'rotation', -30);
        set(gca, 'xcolor', 'k', 'ycolor', 'k');
    else
    text(d+0.1, -0.3, sprintf('p = %.3f', pval), ...
        'fontangle', 'italic', 'fontsize', 4, 'rotation', -30);
        set(gca, 'xcolor', 'k', 'ycolor', 'k');
    end
end

set(gca, 'xtick', 1:length(datasets), 'xticklabel', legtext, 'xticklabelrotation', -30, 'xcolor', 'k');
ylabel({'History shift in v_{bias}'}, 'color', colors(2, :));
axis tight; offsetAxes;
tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/bias_posteriors_driftbias.pdf'));

% ========================================== %
% STARTING POINT POSTERIORS - error vs. correct
% ========================================== %

close; subplot(3,3,1); hold on;
plot([1 6], [0 0], 'k-', 'linewidth', 0.5);

for d = 1:length(datasets),
    disp(datasets{d});
    dat = readtable(sprintf('%s/%s/stimcoding_dc_z_prevcorrect/group_traces.csv', mypath, datasets{d}));
    difference = (dat.z_trans_1_1_ - dat.z_trans_1__1_) - (dat.z_trans_0_1_ - dat.z_trans_0__1_);
        violinPlot_distribution(d, difference, colors2(d, :));
    legtext{d} = cat(2, datasetnames{d}{1}, ' ', datasetnames{d}{2});
    pval = posteriorpval((dat.z_trans_1_1_ - dat.z_trans_1__1_), (dat.z_trans_0_1_ - dat.z_trans_0__1_));
    text(d+0.1, -0.3, sprintf('p = %.3f', pval), ...
        'fontangle', 'italic', 'fontsize', 4, 'rotation', -30);
    set(gca, 'xcolor', 'k', 'ycolor', 'k');
end

set(gca, 'xtick', 1:length(datasets), 'xticklabel', legtext, 'xticklabelrotation', -30, 'xcolor', 'k');
ylabel({'History shift in z'; 'after correct - error'});
axis tight; offsetAxes;
tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/bias_posteriors_startingpoint_outcome.pdf'));


% ========================================== %
% DRIFT BIAS POSTERIORS - error vs. correct
% ========================================== %

close; subplot(3,3,1); hold on;
plot([1 6], [0 0], 'k-', 'linewidth', 0.5);

for d = 1:length(datasets),
    dat = readtable(sprintf('%s/%s/stimcoding_dc_z_prevcorrect/group_traces.csv', mypath, datasets{d}));
    disp(datasets{d});
    difference = (dat.dc_1_1_ - dat.dc_1__1_) -  (dat.dc_0_1_ - dat.dc_0__1_) ;
    violinPlot_distribution(d, difference, colors2(d, :));
    legtext{d} = cat(2, datasetnames{d}{1}, ' ', datasetnames{d}{2});
    pval = posteriorpval((dat.dc_1_1_ - dat.dc_1__1_), (dat.dc_0_1_ - dat.dc_0__1_));
    text(d+0.1, -1, sprintf('p = %.3f', pval), ...
        'fontangle', 'italic', 'fontsize', 4, 'rotation', -30);
        set(gca, 'xcolor', 'k', 'ycolor', 'k');

end

set(gca, 'xtick', 1:length(datasets), 'xticklabel', legtext, 'xticklabelrotation', -30, 'xcolor', 'k');
ylabel({'History shift in v_{bias}'; 'after correct - error'});
axis tight; offsetAxes;
tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/bias_posteriors_driftbias_outcome.pdf'));

end