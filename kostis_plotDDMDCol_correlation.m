function alldat = kostis_plotDDMDCol_correlation

% Code to fit the history-dependent drift diffusion models described in
% Urai AE, Gee JW de, Donner TH (2018) Choice history biases subsequent evidence accumulation. bioRxiv:251595
%
% MIT License
% Copyright (c) Anne Urai, 2018
% anne.urai@gmail.com

close all; clc;
global mypath datasets colors
d = 4;
results = readtable(sprintf('%s/summary/%s/allindividualresults_kostis.csv', mypath, 'Anke_MEG_transition'));

% ============================================ %
% ONE LARGE PLOT WITH PANEL FOR EACH DATASET
% ============================================ %
% PLOT ONLY BAR GRAPH
[rho1, pval1] = corr(results.ddmDColl_dcz_zbias, results.repetitionK, ...
    'type', 'spearman', 'rows', 'complete');

[rho2, pval2] = corr(results.ddmDColl_dcz_dcbias, results.repetitionK, ...
    'type', 'spearman', 'rows', 'complete');

% compute the difference in correlation
[rho3, pval3] = corr(results.ddmDColl_dcz_zbias, results.ddmDColl_dcz_dcbias, ...
    'rows', 'complete', 'type', 'spearman');
[rhodiff, ~, pval] = rddiffci(rho1,rho2,rho3, length(results.repetitionK), 0.05);

subplot(4,6,1);
bar([rho1 rho2], 'facecolor', mean(colors([1 2], :)), 'edgecolor', 'none', 'barwidth', 0.6);
set(gca, 'xtick', 1:2, 'xticklabel', {'z', 'v_{bias}'});
%ylabel({'Correlation \rho', 'with P(repeat'});
mysigstar(gca, 1, 0.05, pval1, 'w');
mysigstar(gca, 2, 0.05, pval2, 'w');
%mysigstar(gca, [1 2],0.95, pval);
axis tight; box off; 
ylim([0 1]);
offsetAxes;
tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/DDMDColl_correlations.pdf'));


end
