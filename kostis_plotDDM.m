function kostis_plotDDM

global mypath datasets datasetnames colors
    
% ============================================ %
% ONE LARGE PLOT WITH PANEL FOR EACH DATASET
% ============================================ %

results = readtable(sprintf('%s/summary/%s/allindividualresults_kostis.csv', mypath, 'Anke_MEG_Neutral'));

% assign to structure
allresults.z            = results.ddmK_dcz_zbias;
allresults.dc           = results.ddmK_dcz_dcbias;
allresults.criterionshift = results.repetition_alldata;

allresults.marker 			= 'o';
allresults.meancolor 		= [0 0 0];
allresults.scattercolor	 	= [0.5 0.5 0.5]
close all;

% PLOT
sp1 = subplot(4,4,1); hold on;
[rho1, tt1] = plotScatter(allresults, 'z', 0.1, 1);
ylabel('P(repeat)');

sp2 = subplot(4,4,2); hold on;
[rho2, tt2, handles] = plotScatter(allresults, 'dc', 0.1, 1);
set(gca, 'yticklabel', []);

set(sp2, 'ylim', get(sp1, 'ylim'), 'ytick', get(sp1, 'ytick'));

% compute the difference in correlation
[rho3, pval3] = corr(cat(1, allresults.z), cat(1, allresults.dc), ...
    'rows', 'complete', 'type', 'pearson');
[rhodiff1, ~, pvalD1] = rddiffci(rho1,rho2,rho3,numel(~isnan( cat(1, allresults(:).criterionshift))), 0.05);
sp2.Position(1) = sp2.Position(1) - 0.08;

% add colored axes after suplabel (which makes them black)

xlabel(sp1, {'History shift in z'});
set(sp1, 'xcolor', colors(1, :), 'ycolor', 'k');

xlabel(sp2, {'History shift in v_{bias}'});
set(sp2, 'xcolor', colors(2, :),  'ycolor', 'k');

%% add line between the two correlation coefficients
txt = {sprintf('\\Deltar(%d) = %.3f, p = %.3f', length(find(~isnan(cat(1, allresults(:).criterionshift) )))-3, rhodiff1, pvalD1)};
if pvalD1 < 0.001,
    txt = {sprintf('\\Deltar(%d) = %.3f, p < 0.001', length(find(~isnan(cat(1, allresults(:).criterionshift) )))-3,  rhodiff1)};
end

tt = title(sp2, txt, 'fontweight', 'normal', 'fontsize', 6, 'horizontalalignment', 'right');

tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/DDM_kostis_correlations.pdf'));

end