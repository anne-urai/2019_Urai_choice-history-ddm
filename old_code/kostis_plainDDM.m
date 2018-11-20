function kostis_plainDDM

global mypath datasetnames colors datasets
d = 2; close all;
results     = readtable(sprintf('%s/summary/%s/allindividualresults_kostis.csv', mypath, 'Anke_MEG_transition'));

% assign to structure
allresults(1).z_prevresp     = results.ddmK_dcz_zbias;
allresults(1).v_prevresp     = results.ddmK_dcz_dcbias;
allresults(1).criterionshift = results.repetitionK;

allresults(1).marker 			= 'o';
allresults(1).meancolor 		= [0 0 0];
allresults(1).scattercolor	 	= [0.5 0.5 0.5];
alltitles{1} 					= {datasetnames{d}{1} datasetnames{d}{2}}; % use only the dataset title

% PLOT
sp1 = subplot(4,4,1); hold on;
[rho1, tt1] = plotScatter(allresults, 'z_prevresp', 0.585, 1);
ylabel('P(repeat)');

sp2 = subplot(4,4,2); hold on;
[rho2, tt2, handles] = plotScatter(allresults, 'v_prevresp', 0.05, 1);
set(gca, 'yticklabel', []);

set(sp2, 'ylim', get(sp1, 'ylim'), 'ytick', get(sp1, 'ytick'));

% compute the difference in correlation
[rho3, pval3] = corr(cat(1, allresults(:).v_prevresp), cat(1, allresults(:).z_prevresp), ...
    'rows', 'complete', 'type', 'spearman');
if pval3 < 0.05,
    fprintf('warning %s: rho = %.3f, pval = %.3f \n', datasets{d}, rho3, pval3);
end
[rhodiff, ~, pval] = rddiffci(rho1,rho2,rho3,numel(~isnan( cat(1, allresults(:).criterionshift))), 0.05);

% move together
sp2.Position(1) = sp2.Position(1) - 0.08;
try
    ss = suplabel(' ', 't');
catch
   % ss = suplabel(datasetnames{d}{1}, 't');
end

% set(ss, 'fontweight', 'normal');
% ss.FontWeight = 'normal';
% ss.Position(2) = ss.Position(2) - 0.03;

% add colored axes after suplabel (which makes them black)
xlabel(sp1, 'History shift in z');
set(sp1, 'xcolor', colors(1, :), 'ycolor', 'k');
xlabel(sp2, 'History shift in v_{bias}');
set(sp2, 'xcolor', colors(2, :), 'ycolor', 'k');

if 1,
    %% add line between the two correlation coefficients
    txt = {sprintf('\\Delta\\rho(%d) = %.3f, p = %.3f', length(find(~isnan(cat(1, allresults(:).criterionshift) )))-3, rhodiff, pval)};
    if pval < 0.001,
        txt = {sprintf('\\Delta\\rho(%d) = %.3f, p < 0.001', length(find(~isnan(cat(1, allresults(:).criterionshift) )))-3,  rhodiff)};
    end
    tt = title(txt, 'fontweight', 'normal', 'fontsize', 6, 'horizontalalignment', 'left');
    tt.Position(2) = tt.Position(2) - 0.008;
end

tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/DDMkostis_plain.pdf'));
end
