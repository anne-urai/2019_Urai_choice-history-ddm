function e11_serialBias_MEGregression_posteriors
% ========================================== %
% conditional response functions from White & Poldrack
% run on simulated rather than real data
% ========================================== %

addpath(genpath('~/code/Tools'));
warning off; close all; clear;
global mypath

datasets        = {'MEG_MEGdata'};
datasetnames    = {'MEG trials'};
d = 1;
mdls = {'regress_dc_z_prevresp', ...
    'regress_dc_z_motorslope', ...
    'regress_dc_z_motorstart', ...
    'regress_dc_z_visualgamma'};

close all;
dat = readtable(sprintf('%s/%s/%s/group_traces.csv', mypath, datasets{d}, 'regress_dc_z_visualgamma'));
subplot(331); plotHist(dat, 'v_visualgamma');
subplot(332); plotHist(dat, 'z_visualgamma');

dat = readtable(sprintf('%s/%s/%s/group_traces.csv', mypath, datasets{d}, 'regress_dc_z_motorslope'));
subplot(334); plotHist(dat, 'v_motorslope');
subplot(335); plotHist(dat, 'z_motorslope');

dat = readtable(sprintf('%s/%s/%s/group_traces.csv', mypath, datasets{d}, 'regress_dc_z_motorstart'));
subplot(337); plotHist(dat, 'v_motorbeta');
subplot(338); plotHist(dat, 'z_motorbeta');

suplabel('Posterior probability', 'y');
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/MEG_posteriors.pdf'));

end

function plotHist(dat, param)
hold on; histogram(dat.(param), 'edgecolor', 'none', 'facecolor', [0.5 0.5 0.5]); vline(0, 'color', 'r');
title(sprintf('p = %.4f', min([mean(dat.(param) < 0) mean(dat.(param) > 0)])));
xlabel(regexprep(param, '_', '~'));
end