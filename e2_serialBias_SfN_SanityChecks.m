addpath(genpath('~/code/Tools'));
warning off; close all; clear;

global datasets datasetnames
usr = getenv('USER');
switch usr
case 'anne' % local
  datasets = {'RT_RDK', 'projects/0/neurodec/Data/MEG-PL', 'NatComm', 'Anke_2afc_neutral'};
case 'aeurai' % lisa/cartesius
  datasets = {'RT_RDK', 'MEG', 'NatComm', 'Anke_neutral', 'Anke_repetitive', 'Anke_alternating'};
end
datasetnames = {'RT', '2IFC',  'NatComm', 'Anke neutral', 'Anke repetitive', 'Anke alternating'};

set(groot, 'defaultaxesfontsize', 6, 'defaultaxestitlefontsizemultiplier', 1.1, ...
    'defaultaxestitlefontweight', 'bold', ...
    'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters', 'defaultaxestickdir', 'out');

% ========================================== %
% DOES DRIFT RATE CORRELATE WITH D'?
% ========================================== %

disp(datasets);

cnt = 1;
for d = 1:length(datasets),
    % cnt = 1+(d-1)*4;

    dat = readtable(sprintf('~/Data/%s/HDDM/summary/allindividualresults.csv', ...
        datasets{d}));
    subplot(4,4,cnt); cnt = cnt + 2;

    if d < 3,

        % find the v_stimulus variable
        varnames = find(~cellfun(@isempty, ...
            regexp(dat.Properties.VariableNames, 'v_stimulus__regressdcprev')));
        % find the one with the most values (i.e. one per session)
        tmpdat          = dat{:, varnames};
        nrdatapoints    = sum(~isnan(tmpdat), 1);
        [~, varidx]     = max(nrdatapoints);
        varname         = dat.Properties.VariableNames{varnames(varidx)};

        % scatterplot
        s = gscatter(dat.dprime, dat.(varname), dat.session, viridis(numel(unique(dat.session))), '.', 5, 'off');
        box off;
        xlabel('d'''); ylabel('v ~ stimulus');
        title(datasetnames{d});

        % layout
        axis tight; axis square;
        xlim([floor(min(get(gca, 'xlim'))) ceil(max(get(gca, 'xlim')))]);
        ylim([floor(min(get(gca, 'ylim'))) ceil(max(get(gca, 'ylim')))]);
        set(gca, 'xtick', min(get(gca, 'xlim')):max(get(gca, 'xlim')));
        set(gca, 'ytick', min(get(gca, 'ylim')):max(get(gca, 'ylim')));

        l = legend(s, {'', 'Session 1', 'Session 2', 'Session 3', 'Session 4', 'Session 5'});
        l.Position(1) = l.Position(1) + 0.12;
        l.Box = 'off';

    else

        % ========================================== %
        % Anke's data - drift rate as a function of coherence
        % ========================================== %

        vars = dat.Properties.VariableNames';
        cohvars = vars(~cellfun(@isempty, strfind(vars, 'dprime_c')));
        alldprime = dat{dat.session == 0, cohvars};
        driftvars = regexp(vars, 'v_c\w+__stimcodingnohist', 'match');
        driftvars = vars((~cellfun(@isempty, driftvars)));
        alldrift  = dat{dat.session == 0, driftvars};

        % which coherences were used?
        cohvars = regexprep(cohvars, '_', '.');
        cohs    = cellfun(@sscanf, cohvars, repmat({'dprime.c%f'}, length(cohvars), 1));
        allcohs = repmat(cohs', size(alldrift, 1), 1);
        g       = gscatter(alldprime(:), alldrift(:), allcohs(:), plasma(numel(unique(allcohs))), '.', 5, 'off');
        xlabel('d'''); ylabel('v ~ stimulus');

        box off;
        l = legend(g, arrayfun(@num2str, cohs, 'un', 0));
        l.Position(1) = l.Position(1) + 0.1;
        l.Box = 'off';
        text(10, 4.5, '% coherence', 'fontsize', 6);

        axis tight; axis square;
        xlim([floor(min(get(gca, 'xlim'))) ceil(max(get(gca, 'xlim')))]);
        ylim([floor(min(get(gca, 'ylim'))) ceil(max(get(gca, 'ylim')))]);
        set(gca, 'xtick', min(get(gca, 'xlim')):max(get(gca, 'xlim')));
        set(gca, 'ytick', min(get(gca, 'ylim')):max(get(gca, 'ylim')));
        title(datasetnames{d});
    end
end
print(gcf, '-depsc', '~/Data/serialHDDM/fig2_driftrate.eps');

% ========================================== %
% DOES DC ~ PREVRESP CORRELATE WITH P(REPEAT)?
% ========================================== %

vars = {'v', 'z'};
for v = 1:length(vars),
    close;
    cnt = 1;
    for d = 1:length(datasets),
        % cnt = 1+(d-1)*4;

        dat = readtable(sprintf('~/Data/%s/HDDM/summary/allindividualresults.csv', ...
            datasets{d}));
        subplot(4,4,cnt); cnt = cnt + 1;

        % find the v_stimulus variable
        varnames = find(~cellfun(@isempty, ...
            regexp(dat.Properties.VariableNames, sprintf('%s_prevresp__regressdczprev', vars{v}))));

        % find the one with the most values (i.e. one per session)
        tmpdat          = dat{:, varnames};
        nrdatapoints    = sum(~isnan(tmpdat), 1);
        [~, varidx]     = max(nrdatapoints);
        varname         = dat.Properties.VariableNames{varnames(varidx)};

        % switch sign for z!
        switch vars{v}
            case 'z'
                dat.(varname) = -1 * dat.(varname);
        end

        % scatterplot
        s = scatter(dat.repetition, dat.(varname), 50, dat.session, '.');
        box off;
        xlabel('p(repeat)'); ylabel(sprintf('%s ~ prevresp', vars{v}));
        title(datasetnames{d});

        [rho, pval] = corr(dat.repetition, dat.(varname), 'type', 'spearman', 'rows', 'complete');
        if pval < 0.05, lsline; end;

        axis square;
        xlim([0 1]); ylim([-1 1]);

        text(0.1, -0.6, ...
            {sprintf('\\rho = %.2f', rho); sprintf('p = %.3f', pval)}, 'fontsize', 5);

        if d == 3, cnt = cnt + 1; end
    end
    print(gcf, '-depsc', sprintf('~/Data/serialHDDM/fig2_repetitionprobability_%s.eps', vars{v}));
end

% ========================================== %
% MODULATION OF BOTH DC AND Z
% ========================================== %

close;
cnt = 1;
for d = 1:length(datasets),

    dat = readtable(sprintf('~/Data/%s/HDDM/summary/allindividualresults.csv', ...
        datasets{d}));
    subplot(4,4,cnt); cnt = cnt + 1;

    for v = 1:length(vars),
        % find the v_stimulus variable
        varnames = find(~cellfun(@isempty, ...
            regexp(dat.Properties.VariableNames, sprintf('%s_prevresp__regressdczprev', vars{v}))));

        % find the one with the most values (i.e. one per session)
        tmpdat          = dat{:, varnames};
        nrdatapoints    = sum(~isnan(tmpdat), 1);
        [~, varidx]     = max(nrdatapoints);
        varname_both{v}         = dat.Properties.VariableNames{varnames(varidx)};
    end
    % scatterplot

    s = scatter(dat.(varname_both{1}), dat.(varname_both{2}), 50, dat.session, '.');
    box off;
    xlabel(sprintf('%s ~ prevresp', vars{1})); ylabel(sprintf('%s ~ prevresp', vars{2}));
    title(datasetnames{d});

    [rho, pval] = corr(dat.(varname_both{1}), dat.(varname_both{2}), 'type', 'spearman', 'rows', 'complete');
    if pval < 0.05, lsline; end;

    axis square;
  %  text(0.1, -0.6, ...
    %    {sprintf('\\rho = %.2f', rho); sprintf('p = %.3f', pval)}, 'fontsize', 5);

    if d == 3, cnt = cnt + 1; end
end
print(gcf, '-depsc', sprintf('~/Data/serialHDDM/fig2_repetitionprobability_dcVz.eps'));

% ========================================== %
% TEST CORRELATION BETWEEN DC AND Z CORR
% ========================================== %

switch usr
    case 'anne' % local
        datasets = { 'Anke_2afc_neutral', 'Anke_2afc_repetitive', 'Anke_2afc_alternating'};
    case 'aeurai' % lisa/cartesius
        datasets = { 'Anke_neutral', 'Anke_repetitive', 'Anke_alternating'};
end

for d = 1:length(datasets),
    dat = readtable(sprintf('~/Data/%s/HDDM/summary/allindividualresults.csv', ...
        datasets{d}));

    var1 = dat.v_prevresp__regressdczprevrespstim;
    var2 = -dat.z_prevresp__regressdczprevrespstim;
    var3 = dat.repetition;
    r12 = corr(var1, var2, 'rows', 'complete');
    r13 = corr(var1, var3, 'rows', 'complete');
    r23 = corr(var2, var3, 'rows', 'complete');

    [rhodiff(d),ci, pval(d)] = rddiffci(r13, r23, r12, 22, 0.05);

end
