function alldat = kostis_plotRamp_correlation

% Code to fit the history-dependent drift diffusion models described in
% Urai AE, Gee JW de, Donner TH (2018) Choice history biases subsequent evidence accumulation. bioRxiv:251595
%
% MIT License
% Copyright (c) Anne Urai, 2018
% anne.urai@gmail.com

close all; clc;
global mypath datasets
results = readtable(sprintf('%s/summary/%s/allindividualresults_kostis.csv', mypath, 'Anke_MEG_transition'));

cols1 = cbrewer('qual', 'Set1', 8);
cols2 = cbrewer('qual', 'Dark2', 8);
thiscol = nanmean([cols1(2, :); cols2(6, :)]);

% ============================================ %
% ONE LARGE PLOT WITH PANEL FOR EACH DATASET
% ============================================ %

% PLOT ONLY BAR GRAPH
[rho1, pval1] = corr(results.ddmK_rp2_offset, results.repetitionK, ...
    'type', 'spearman', 'rows', 'complete');

[rho2, pval2] = corr(results.ddmK_rp2_slope, results.repetitionK, ...
    'type', 'spearman', 'rows', 'complete');

% compute the difference in correlation
[rho3, pval3] = corr(results.ddmK_rp2_offset, results.ddmK_rp2_slope, ...
    'rows', 'complete', 'type', 'spearman');
[rhodiff, ~, pval] = rddiffci(rho1,rho2,rho3, length(results.repetitionK), 0.05);

subplot(4,6,1);
bar([rho1 rho2], 'facecolor', thiscol, 'edgecolor', 'none', 'barwidth', 0.6);
set(gca, 'xtick', 1:2, 'xticklabel', {'offset', 'ramp'});
% ylabel({'Correlation \rho', 'with P(repeat'});
%ylabel({'Correlation \rho', 'with P(repeat'});
mysigstar(gca, 1, 0.05, pval1, 'w');
mysigstar(gca, 2, 0.05, pval2, 'w');
mysigstar(gca, [1 2],0.95, pval);
axis tight; box off; 
ylim([0 1]);
offsetAxes;
tightfig;

print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/DDMramp_correlations.pdf'));


end
