function alldat = kostis_plotRamp_correlation

% Code to fit the history-dependent drift diffusion models described in
% Urai AE, Gee JW de, Donner TH (2018) Choice history biases subsequent evidence accumulation. bioRxiv:251595
%
% MIT License
% Copyright (c) Anne Urai, 2018
% anne.urai@gmail.com

close all; clc;
global mypath datasets datasetnames colors

results = readtable(sprintf('%s/summary/%s/allindividualresults_kostis.csv', mypath, 'Anke_MEG_Neutral'));

% ============================================ %
% ONE LARGE PLOT WITH PANEL FOR EACH DATASET
% ============================================ %

results.criterionshift = results.repetition_alldata;

% assign to structure
allresults.z_prevresp     = results.ouK_sp_spbias;
allresults.v_prevresp     = results.ouK_input_inputbias;
allresults.criterionshift = results.criterionshift;

allresults.marker 			= 'o';
allresults.meancolor 		= [ 0 0 0];
allresults.scattercolor	 	= [ 0.5588    0.5824    0.3647];
close all;

% PLOT
sp1 = subplot(4,4,1); hold on;
[rho1, tt1] = plotScatter(allresults, 'z_prevresp', 0.1, 1);
ylabel('P(repeat)');

sp2 = subplot(4,4,2); hold on;
[rho2, tt2, handles] = plotScatter(allresults, 'v_prevresp', 0.6, 1);
set(gca, 'yticklabel', []);

set(sp2, 'ylim', get(sp1, 'ylim'), 'ytick', get(sp1, 'ytick'));

% compute the difference in correlation
[rho3, pval3] = corr(cat(1, allresults(:).v_prevresp), cat(1, allresults(:).z_prevresp), ...
    'rows', 'complete', 'type', 'pearson');
if pval3 < 0.05,
    fprintf('warning %s: rho = %.3f, pval = %.3f \n', datasets{d}, rho3, pval3);
end
[rhodiff, ~, pval] = rddiffci(rho1,rho2,rho3,numel(~isnan( cat(1, allresults(:).criterionshift))), 0.05);

% move together
% set(ss, 'fontweight', 'normal');
% ss.FontWeight = 'normal';
% ss.Position(2) = ss.Position(2) - 0.03;
sp2.Position(1) = sp2.Position(1) - 0.08;

% add colored axes after suplabel (which makes them black)

xlabel(sp2, {'Ramp slope'});
set(sp2, 'xcolor', 'k',  'ycolor', 'k');
xlabel(sp1, {'Ramp offset'});
set(sp1, 'xcolor', 'k', 'ycolor', 'k');
%% add line between the two correlation coefficients
txt = {sprintf('\\Deltar(%d) = %.3f, p = %.3f', length(find(~isnan(cat(1, allresults(:).criterionshift) )))-3, rhodiff, pval)};
if pval < 0.001,
    txt = {sprintf('\\Deltar(%d) = %.3f, p < 0.001', length(find(~isnan(cat(1, allresults(:).criterionshift) )))-3,  rhodiff)};
end

tt = title(sp1, txt, 'fontweight', 'normal', 'fontsize', 6, 'horizontalalignment', 'left');
% tt.Position(2) = tt.Position(2) - 0.008;

tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/DDMramp_correlations.pdf'));


end
