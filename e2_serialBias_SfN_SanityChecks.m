addpath(genpath('~/code/Tools'));
warning off; close all; clear;
global datasets datasetnames
datasets = {'RT_RDK', 'projects/0/neurodec/Data/MEG-PL', 'Anke_2afc_neutral', 'Anke_2afc_sequential'};
datasetnames = {'RT', '2IFC', 'Anke neutral', 'Anke all'};

set(groot, 'defaultaxesfontsize', 6, 'defaultaxestitlefontsizemultiplier', 1.1, ...
    'defaultaxestitlefontweight', 'bold', ...
    'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters');

% ========================================== %
% DOES DRIFT RATE CORRELATE WITH D'?
% ========================================== %

cnt = 1;
for d = 1:2,
    
    dat = readtable(sprintf('~/Data/%s/HDDM/summary/allindividualresults.csv', ...
        datasets{d}));
    subplot(4,4,cnt); cnt = cnt + 1;
    
    % find the v_stimulus variable
    varnames = find(~cellfun(@isempty, ...
        regexp(dat.Properties.VariableNames, 'v_stimulus__regressdcprev')));
    
    % find the one with the most values (i.e. one per session)
    tmpdat          = dat{:, varnames};
    nrdatapoints    = sum(~isnan(tmpdat), 1);
    [~, varidx]     = max(nrdatapoints);
    varname         = dat.Properties.VariableNames{varnames(varidx)};
    
    % scatterplot
    s = gscatter(dat.dprime, dat.(varname), dat.session, viridis(numel(unique(dat.session))), '.', 8, 'off');
    box off;
    xlabel('d'''); ylabel('v ~ stimulus');
    title(datasetnames{d});
    
    % layout
    axis tight; axis square;
    xlim([floor(min(get(gca, 'xlim'))) ceil(max(get(gca, 'xlim')))]);
    ylim([floor(min(get(gca, 'ylim'))) ceil(max(get(gca, 'ylim')))]);
    set(gca, 'xtick', min(get(gca, 'xlim')):max(get(gca, 'xlim')));
    set(gca, 'ytick', min(get(gca, 'ylim')):max(get(gca, 'ylim')));
    
    % show that dc_prevresp correlates with p_repeat
    subplot(4,4,cnt); cnt = cnt + 1;
    
    % find the v_stimulus variable
    varnames = find(~cellfun(@isempty, ...
        regexp(dat.Properties.VariableNames, 'v_prevresp__regressdcprev')));
    
    % find the one with the most values (i.e. one per session)
    tmpdat          = dat{:, varnames};
    nrdatapoints    = sum(~isnan(tmpdat), 1);
    [~, varidx]     = max(nrdatapoints);
    varname         = dat.Properties.VariableNames{varnames(varidx)};
    
    % scatterplot
    s = gscatter(dat.repetition, dat.(varname), dat.session, viridis(numel(unique(dat.session))), '.', 8, 'off');
    box off;
    xlabel('p(repeat)'); ylabel('v ~ prevresp');
    title(datasetnames{d});
    
    % layout
    %axis tight;
    axis square;
    xlim([0.4 0.6]); ylim([-0.3 0.3]);
    set(gca, 'ytick', [-0.3:0.1:0.3]);
    
    l = legend(s, {'', 'Session 1', 'Session 2', 'Session 3', 'Session 4', 'Session 5'});
    l.Position(1) = l.Position(1) + 0.15;
    l.Box = 'off';
    
    cnt = cnt + 2;
end

% ========================================== %
% Anke's data - drift rate as a function of coherence
% ========================================== %

d = 3;
dat = readtable(sprintf('~/Data/%s/HDDM/summary/allindividualresults.csv', ...
    datasets{d}));
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
g = gscatter(alldprime(:), alldrift(:), allcohs(:), plasma(numel(unique(allcohs))), '.', 8, 'off');
xlabel('d'''); ylabel('v ~ stimulus');


box off;
l = legend(g, {'0', '5', '10', '20', '40', '60'});
l.Position(1) = l.Position(1) + 0.08;
l.Box = 'off';
text(6, 4.5, '% coherence', 'fontsize', 6);

axis tight; axis square;
xlim([floor(min(get(gca, 'xlim'))) ceil(max(get(gca, 'xlim')))]);
ylim([floor(min(get(gca, 'ylim'))) ceil(max(get(gca, 'ylim')))]);
set(gca, 'xtick', min(get(gca, 'xlim')):max(get(gca, 'xlim')));
set(gca, 'ytick', min(get(gca, 'ylim')):max(get(gca, 'ylim')));
title(datasetnames{d});

% ======================================================= %
% ANKE'S DATA - correlate accuracy with serial bias
% ======================================================= %

d = 4;
dat = readtable(sprintf('~/Data/%s/HDDM/summary/allindividualresults.csv', ...
    datasets{d}));

% get the sessions separately per condition
subplot(4,4,11);
allaccuracy     = [splitapply(@nanmean, dat.accuracy(dat.transitionprob == 0.2), findgroups(dat.subjnr(dat.transitionprob == 0.2))),  ...
    splitapply(@nanmean, dat.accuracy(dat.transitionprob == 0.5), findgroups(dat.subjnr(dat.transitionprob == 0.5))) ...
    splitapply(@nanmean, dat.accuracy(dat.transitionprob == 0.8), findgroups(dat.subjnr(dat.transitionprob == 0.8)))];

allserialbias   = dat{dat.session == 0, {'v_prevresp_alternating__regressdcprevresp', ...
    'v_prevresp_neutral__regressdcprevresp', ...
    'v_prevresp_repetitive__regressdcprevresp'}};

allcondition = repmat([0.2 0.5 0.8], size(allserialbias, 1), 1);
g = gscatter(allserialbias(:), allaccuracy(:), allcondition(:), flipud(linspecer(3)), '.', 8, 'off');
xlabel('v ~ prevresp'); ylabel('accuracy (%)');
lsline; axisNotSoTight; axis square; box off;
xlim([-1 1]); ylim([0.6 1]); set(gca, 'ytick', 0.5:0.1:1);

% make correlations
for i = 1:3,
    [rho(i), pval(i)] = corr(allaccuracy(:, i), allserialbias(:, i), 'type', 'spearman');
end
l = legend(g, {sprintf('Repetitive, \\rho = %.2f, p = %.3f', rho(1), pval(1)), ...
    sprintf('Neutral, \\rho = %.2f, p = %.3f', rho(2), pval(2)),...
    sprintf('Alternating, \\rho = %.2f, p = %.3f', rho(3), pval(3))});
l.Position(1) = l.Position(1) + 0.25;
l.Box = 'off';
title(datasetnames{d});

print(gcf, '-dpdf', '~/Data/serialHDDM/fig2_driftrate.pdf');


