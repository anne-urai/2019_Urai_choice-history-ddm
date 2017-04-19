addpath(genpath('~/code/Tools'));
warning off; close all; clear;
global datasets datasetnames
datasets = {'RT_RDK', 'projects/0/neurodec/Data/MEG-PL', 'Anke_2afc_sequential'};
datasetnames = {'RT', '2IFC', 'Anke'};

set(groot, 'defaultaxesfontsize', 7, 'defaultaxestitlefontsizemultiplier', 1.1, ...
    'defaultaxestitlefontweight', 'bold', ...
    'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters');

% ========================================== %
% DOES DRIFT RATE CORRELATE WITH D'?
% ========================================== %

for d = 1:length(datasets),
    
    dat = readtable(sprintf('~/Data/%s/HDDM/summary/allindividualresults.csv', ...
        datasets{d}));
    subplot(4,4,d);
    
    % find the v_stimulus variable
    varnames = find(~cellfun(@isempty, ...
        regexp(dat.Properties.VariableNames, 'v_stimulus__regressdcprev')));
    
    % find the one with the most values (i.e. one per session)
    tmpdat          = dat{:, varnames};
    nrdatapoints    = sum(~isnan(tmpdat), 1);
    [~, varidx]     = max(nrdatapoints);
    varname         = dat.Properties.VariableNames{varnames(varidx)};
    
    % scatterplot
    scatter(dat.dprime, ...
        dat.(varname), ...
        8, dat.session);
    
    xlabel('d'''); ylabel('v');
    title(datasetnames{d});
    
    % layout
    axis tight; axis square; lsline;
    xlim([floor(min(get(gca, 'xlim'))) ceil(max(get(gca, 'xlim')))]);
    ylim([floor(min(get(gca, 'ylim'))) ceil(max(get(gca, 'ylim')))]);
    set(gca, 'xtick', min(get(gca, 'xlim')):max(get(gca, 'xlim')));
    set(gca, 'ytick', min(get(gca, 'ylim')):max(get(gca, 'ylim')));
    
end

% ========================================== %
% Anke's data - drift rate as a function of coherence
% ========================================== %

subplot(4,4,9);
alldprime = dat{dat.session == 0, {'dprime_c0', 'dprime_c5', 'dprime_c10', ...
    'dprime_c20', 'dprime_c40', 'dprime_c60'}};
alldrift = dat{dat.session == 0, {'v_c0__stimcodingdcprevrespstim', ...
    'v_c5__stimcodingdcprevrespstim', ...
    'v_c10__stimcodingdcprevrespstim', ...
    'v_c20__stimcodingdcprevrespstim', ...
    'v_c40__stimcodingdcprevrespstim', ...
    'v_c60__stimcodingdcprevrespstim'}};
allcohs = repmat([0 5 10 20 40 60], size(alldrift, 1), 1);
g = gscatter(alldprime(:), alldrift(:), allcohs(:), viridis(6), '.', 8, 'off');
xlabel('d'''); ylabel('v');

box off;
l = legend(g, {'0', '5', '10', '20', '40', '60'});
l.Position(1) = l.Position(1) + 0.08;
l.Box = 'off';
text(6, 4.5, '% coherence', 'fontsize', 7);

axis tight; axis square;
xlim([floor(min(get(gca, 'xlim'))) ceil(max(get(gca, 'xlim')))]);
ylim([floor(min(get(gca, 'ylim'))) ceil(max(get(gca, 'ylim')))]);
set(gca, 'xtick', min(get(gca, 'xlim')):max(get(gca, 'xlim')));
set(gca, 'ytick', min(get(gca, 'ylim')):max(get(gca, 'ylim')));
title(datasetnames{d});

% ======================================================= %
% ANKE'S DATA - correlate accuracy with serial bias
% ======================================================= %

% get the sessions separately per condition
subplot(4,4,11);
allaccuracy     = [splitapply(@nanmean, dat.accuracy(dat.transitionprob == 0.2), findgroups(dat.subjnr(dat.transitionprob == 0.2))),  ...
    splitapply(@nanmean, dat.accuracy(dat.transitionprob == 0.5), findgroups(dat.subjnr(dat.transitionprob == 0.5))) ...
    splitapply(@nanmean, dat.accuracy(dat.transitionprob == 0.8), findgroups(dat.subjnr(dat.transitionprob == 0.8)))];

allserialbias   = dat{dat.session == 0, {'v_prevresp_alternating__regressdcprevresp', ...
    'v_prevresp_neutral__regressdcprevresp', ...
    'v_prevresp_repetitive__regressdcprevresp'}};

allcondition = repmat([0.2 0.5 0.8], size(allserialbias, 1), 1);
g = gscatter(allserialbias(:), allaccuracy(:), allcondition(:), (linspecer(3)), '.', 8, 'off');
xlabel('Serial bias (dc)'); ylabel('Accuracy (%)');
lsline; axisNotSoTight; axis square; box off;
xlim([-1 1]); ylim([0.6 1]); set(gca, 'ytick', 0.5:0.1:1);

% make correlations
for i = 1:3,
    [rho(i), pval(i)] = corr(allaccuracy(:, i), allserialbias(:, i), 'type', 'spearman');
end
l = legend(g, {sprintf('alternating, \\rho = %.2f, p = %.3f', rho(1), pval(1)), ...
    sprintf('neutral, \\rho = %.2f, p = %.3f', rho(2), pval(2)),...
    sprintf('repetitive, \\rho = %.2f, p = %.3f', rho(3), pval(3))});
l.Position(1) = l.Position(1) + 0.25;
l.Box = 'off';

print(gcf, '-dpdf', '~/Data/serialHDDM/fig2_driftrate.pdf');


