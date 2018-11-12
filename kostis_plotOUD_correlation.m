function kostis_plotOUD_correlation

global mypath 

%% compare correlation coefficients
close all; 
results = readtable(sprintf('%s/summary/%s/allindividualresults_kostis.csv', mypath, 'Anke_MEG_transition'));
% corrplot(results, {'ouD_input_inputbias', 'ouD_lambda_lambdabias', 'ouD_sp_spbias', 'repetitionK'});
% print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/OU_scatter.pdf'));
% 
% close; 
% corrplot(results, {'repetitionKD', 'ouD_input_inputbias', 'ouD_lambda_lambdabias', 'ouD_sp_spbias'});
% print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/OU_scatter_2.pdf'));

%% separately, test correlation coefficients
    
y = [results.ouD_input_inputbias results.ouD_lambda_lambdabias];

% ============================================ %
% ONE LARGE PLOT WITH PANEL FOR EACH DATASET
% ============================================ %

results.criterionshift = results.repetitionK;

% assign to structure
allresults.sp       = results.ouD_sp_spbias;
allresults.input    = results.ouD_input_inputbias;
allresults.leak     = results.ouD_lambda_lambdabias;
allresults.criterionshift = results.criterionshift;

cols2 = cbrewer('qual', 'Dark2', 8);
colors = cols2([5 3 4], :);

allresults.marker 			= 'o';
allresults.meancolor 		= [ 0 0 0];
allresults.scattercolor	 	= colors(1, :);
close all;

% PLOT
sp1 = subplot(4,4,1); hold on;
[rho1, tt1] = plotScatter(allresults, 'sp', 0.1, 1);
ylabel('P(repeat)');

sp2 = subplot(4,4,2); hold on;
allresults.scattercolor = colors(2, :);
[rho2, tt2, handles] = plotScatter(allresults, 'input', 0.1, 1);
set(gca, 'yticklabel', []);

sp3 = subplot(4,4,3); hold on;
allresults.scattercolor = colors(3, :);
[rho3, tt3, handles] = plotScatter(allresults, 'leak', 0.1, 1);
set(gca, 'yticklabel', []);

set(sp2, 'ylim', get(sp1, 'ylim'), 'ytick', get(sp1, 'ytick'));
set(sp3, 'ylim', get(sp1, 'ylim'), 'ytick', get(sp1, 'ytick'));

% compute the difference in correlation
[rho3, pval3] = corr(cat(1, allresults.sp), cat(1, allresults.input), ...
    'rows', 'complete', 'type', 'pearson');
[rhodiff1, ~, pvalD1] = rddiffci(rho1,rho2,rho3,numel(~isnan( cat(1, allresults(:).criterionshift))), 0.05);


[rho3, pval3] = corr(cat(1, allresults.input), cat(1, allresults.leak), ...
    'rows', 'complete', 'type', 'pearson');
[rhodiff2, ~, pvalD2] = rddiffci(rho1,rho2,rho3,numel(~isnan( cat(1, allresults(:).criterionshift))), 0.05);

% move together
% set(ss, 'fontweight', 'normal');
% ss.FontWeight = 'normal';
% ss.Position(2) = ss.Position(2) - 0.03;
sp2.Position(1) = sp2.Position(1) - 0.08;
sp3.Position(1) = sp3.Position(1) - 0.16;

% add colored axes after suplabel (which makes them black)


xlabel(sp1, {'Offset bias'});
set(sp1, 'xcolor', 'k', 'ycolor', 'k');

xlabel(sp2, {'Input bias'});
set(sp2, 'xcolor', 'k',  'ycolor', 'k');


xlabel(sp3, {'Leak bias'});
set(sp3, 'xcolor', 'k',  'ycolor', 'k');

%% add line between the two correlation coefficients
txt = {sprintf('\\Delta\\rho(%d) = %.3f, p = %.3f', length(find(~isnan(cat(1, allresults(:).criterionshift) )))-3, rhodiff1, pvalD1)};
if pvalD1 < 0.001,
    txt = {sprintf('\\Delta\\rho(%d) = %.3f, p < 0.001', length(find(~isnan(cat(1, allresults(:).criterionshift) )))-3,  rhodiff1)};
end

tt = title(sp1, txt, 'fontweight', 'normal', 'fontsize', 6, 'horizontalalignment', 'center');

txt = {sprintf('\\Delta\\rho(%d) = %.3f, p = %.3f', length(find(~isnan(cat(1, allresults(:).criterionshift) )))-3, rhodiff2, pvalD2)};
if pvalD2 < 0.001,
    txt = {sprintf('\\Delta\\rho(%d) = %.3f, p < 0.001', length(find(~isnan(cat(1, allresults(:).criterionshift) )))-3,  rhodiff2)};
end

tt = title(sp3, txt, 'fontweight', 'normal', 'fontsize', 6, 'horizontalalignment', 'center');

% tt.Position(2) = tt.Position(2) - 0.008;

tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/OUD_correlations.pdf'));


%% =========================================== %
% OVERALL LEAK PARAMETER
% ============================================ %

results = readtable(sprintf('%s/summary/%s/allindividualresults_kostis.csv', mypath, 'Anke_MEG_transition'));

close all;
sp = subplot(3,5,1);
hold on;
% plot identity line
plot([0.8 1.2], [0 0], '-', 'linewidth', 0.5, 'color', [0.5 0.5 0.5]);

lambda = results.ouD_lambda_lambda;
h = scatter(1*ones(size(lambda)), lambda, 3, colors(3, :), 'jitter', 'on', 'jitteramount', 0.07);
plot([1-0.1 1+0.1], [nanmean(lambda) nanmean(lambda)], 'k-');
pval = permtest(lambda);
[h, pval] = ttest(lambda)
mysigstar(gca, 1, max(get(gca, 'ylim')), pval);
%legtext{d} = cat(2, datasetnames{d}{1}, ' ', datasetnames{d}{2});
ylabel('Effective leak (\lambda)');
set(gca, 'xtick', [0.9 1.1], 'xticklabel', []); offsetAxes;
title({'Leaky accumulator', ''});
tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/OUD_lambda.pdf'));

%% ============================================== %
% NOW INCLUDING COLLAPSING BOUND
% ============================================== %
results = readtable(sprintf('%s/summary/%s/allindividualresults_kostis.csv', mypath, 'Anke_MEG_transition'));

close all;
sp = subplot(3,5,1);
hold on;
% plot identity line
plot([0.8 1.2], [0 0], '-', 'linewidth', 0.5, 'color', [0.5 0.5 0.5]);

lambda = results.ouDColl_input_lambda;
h = scatter(1*ones(size(lambda)), lambda, 3, colors(3, :), 'jitter', 'on', 'jitteramount', 0.07);
plot([1-0.1 1+0.1], [nanmean(lambda) nanmean(lambda)], 'k-');
pval = permtest(lambda);
[h, pval] = ttest(lambda)
mysigstar(gca, 1, 0.95*max(get(gca, 'ylim')), pval);
%legtext{d} = cat(2, datasetnames{d}{1}, ' ', datasetnames{d}{2});
%ylabel('Effective leak (\lambda)');
title({'Leaky accumulator', 'collapsing bounds'});

set(gca, 'xtick', [0.9 1.1], 'xticklabel', []); offsetAxes;
tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/OUDColl_lambda.pdf'));

end