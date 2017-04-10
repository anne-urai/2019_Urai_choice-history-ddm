% from the csv table, make an overview of repetition behaviour

% get a huge list with values for each participant
% can then work with this dataframe

clear; close all; clc;
addpath(genpath('~/code/Tools'));
warning off;

datasets = 'Anke_2afc_sequential/HDDM';
set(groot, 'defaultaxesfontsize', 8, 'defaultaxestitlefontsizemultiplier', 1, ...
    'defaultaxestitlefontweight', 'bold', ...
    'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters');
colormap viridis;

results = readtable(sprintf('~/Data/%s/allindividualresults.csv', datasets));
results.Properties.VariableNames{'v_prevresp__regress_dc_z_prevresp'}       = 'dc_seq_regress_joint';
results.Properties.VariableNames{'z_prevresp__regress_dc_z_prevresp'}       = 'z_seq_regress_joint';
results.Properties.VariableNames{'v_prevrt_prevresp__regress_dc_z_prevresp_prevpupil_prevrt'} ...
    = 'dc_rt_seq_regress_joint';
results.Properties.VariableNames{'z_prevrt_prevresp__regress_dc_z_prevresp_prevpupil_prevrt'} ...
    = 'z_rt_seq_regress_joint';


% recode starting point values
z_link_func = @(x) 1 ./ (1 + exp(x));
results.z_seq_regress_joint = ...
    z_link_func(results.z_seq_regress_joint);
results.z_rt_seq_regress_joint = ...
    z_link_func(results.z_rt_seq_regress_joint);

% recode drift criterion into positive = more repetition

% ============================================ %
% BETWEEN THE DIFFERENT TRANSITION PROBABILITY SESSIONS,
% HOW DO PEOPLE BEHAVE?
% ============================================ %

transprob = [0.2 0.5 0.8];

clf;
subplot(331); hold on;
dc = [results.dc_seq_regress_joint(results.transprob == 0.2)'
    results.dc_seq_regress_joint(results.transprob == 0.5)'
    results.dc_seq_regress_joint(results.transprob == 0.8)']';
tp        = repmat(transprob, size(dc, 1), 1);
plot(transprob, dc, 'color', [0.8 0.8 0.8]);
scatter(tp(:), dc(:), 10, tp(:));
xlabel('Transition probability'); ylabel('Drift criterion'); axis square;
box off; axisNotSoTight; set(gca, 'xtick', transprob);
title('Prevresp');

subplot(332); hold on;
z = [results.z_seq_regress_joint(results.transprob == 0.2)'
    results.z_seq_regress_joint(results.transprob == 0.5)'
    results.z_seq_regress_joint(results.transprob == 0.8)']';
plot(transprob, z, 'color', [0.8 0.8 0.8]);
scatter(tp(:), z(:), 10, tp(:));
xlabel('Transition probability'); ylabel('Starting point'); axis square;
box off; axisNotSoTight; set(gca, 'xtick', transprob);
title('Prevresp');

subplot(333);
scatter(results.dc_seq_regress_joint, results.z_seq_regress_joint, ...
    10, results.transprob);
axis square; axisNotSoTight;
xlabel('Drift criterion'); ylabel('Starting point');
lsline; title('Prevresp');
title('Prevresp');

% now the RT addition
subplot(334); hold on;
dc = [results.dc_rt_seq_regress_joint(results.transprob == 0.2)'
    results.dc_rt_seq_regress_joint(results.transprob == 0.5)'
    results.dc_rt_seq_regress_joint(results.transprob == 0.8)']';
plot(transprob, dc, 'color', [0.8 0.8 0.8]);
scatter(tp(:), dc(:), 10, tp(:));
xlabel('Transition probability'); ylabel('Drift criterion'); axis square;
box off; axisNotSoTight; set(gca, 'xtick', transprob);
title('RT x prevresp');

subplot(335); hold on;
z = [results.z_rt_seq_regress_joint(results.transprob == 0.2)'
    results.z_rt_seq_regress_joint(results.transprob == 0.5)'
    results.z_rt_seq_regress_joint(results.transprob == 0.8)']';
plot(transprob, z, 'color', [0.8 0.8 0.8]);
scatter(tp(:), z(:), 10, tp(:));
xlabel('Transition probability'); ylabel('Starting point'); axis square;
box off; axisNotSoTight; set(gca, 'xtick', transprob);
title('RT x prevresp');

subplot(336);
scatter(results.dc_rt_seq_regress_joint, results.z_rt_seq_regress_joint, ...
    10, results.transprob);
axis square; axisNotSoTight;
xlabel('Drift criterion'); ylabel('Starting point');
title('RT x prevresp');

print(gcf, '-dpdf', sprintf('~/Data/%s/Figures/HDDMoverview.pdf', datasets));
