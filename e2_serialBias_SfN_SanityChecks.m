addpath(genpath('~/code/Tools'));
warning off; close all; clear;
global datasets datasetnames
datasets = {'RT_RDK', 'projects/0/neurodec/Data/MEG-PL', 'Anke_2afc_neutral', 'NatComm'};
datasetnames = {'RT', '2IFC', 'Anke neutral', 'NatComm'};

set(groot, 'defaultaxesfontsize', 6, 'defaultaxestitlefontsizemultiplier', 1.1, ...
    'defaultaxestitlefontweight', 'bold', ...
    'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters');

% ========================================== %
% DOES DRIFT RATE CORRELATE WITH D'?
% ========================================== %

cnt = 1;
for d = 1:length(datasets),
    cnt = 1+(d-1)*4;
    
    dat = readtable(sprintf('~/Data/%s/HDDM/summary/allindividualresults.csv', ...
        datasets{d}));
    subplot(4,4,cnt);
    
    if d < 3,
        
        % find the v_stimulus variable
        varnames = find(~cellfun(@isempty, ...
            regexp(dat.Properties.VariableNames, 'v_stimulus__regressdcprevrespstimrt')));
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
        
        l = legend(s, {'', 'Session 1', 'Session 2', 'Session 3', 'Session 4', 'Session 5'});
        l.Position(1) = l.Position(1) + 0.12;
        l.Box = 'off';
        cnt = cnt + 2;
        
    else
        
        % ========================================== %
        % Anke's data - drift rate as a function of coherence
        % ========================================== %
        
        vars = dat.Properties.VariableNames';
        cohvars = vars(~cellfun(@isempty, strfind(vars, 'dprime_c')));
        alldprime = dat{dat.session == 0, cohvars};
        driftvars = regexp(vars, 'v_c\w+__stimcodingdcprevrespstim', 'match');
        driftvars = vars((~cellfun(@isempty, driftvars)));
        alldrift  = dat{dat.session == 0, driftvars};
        
        % which coherences were used?
        cohvars = regexprep(cohvars, '_', '.');
        cohs    = cellfun(@sscanf, cohvars, repmat({'dprime.c%f'}, length(cohvars), 1));
        allcohs = repmat(cohs', size(alldrift, 1), 1);
        g       = gscatter(alldprime(:), alldrift(:), allcohs(:), plasma(numel(unique(allcohs))), '.', 8, 'off');
        xlabel('d'''); ylabel('v ~ stimulus');
        
        box off;
        l = legend(g, arrayfun(@num2str, cohs, 'un', 0));
        l.Position(1) = l.Position(1) + 0.1;
        l.Box = 'off';
        text(7, 4.5, '% coherence', 'fontsize', 6);
        
        axis tight; axis square;
        xlim([floor(min(get(gca, 'xlim'))) ceil(max(get(gca, 'xlim')))]);
        ylim([floor(min(get(gca, 'ylim'))) ceil(max(get(gca, 'ylim')))]);
        set(gca, 'xtick', min(get(gca, 'xlim')):max(get(gca, 'xlim')));
        set(gca, 'ytick', min(get(gca, 'ylim')):max(get(gca, 'ylim')));
        title(datasetnames{d});
        cnt = cnt + 2;
        
    end
    
    % show that dc_prevresp correlates with p_repeat
    subplot(4,4,cnt);
    
    % find the v_stimulus variable
    varnames = find(~cellfun(@isempty, ...
        regexp(dat.Properties.VariableNames, 'v_prevresp__regressdcprev')));
    
    % find the one with the most values (i.e. one per session)
    tmpdat          = dat{:, varnames};
    nrdatapoints    = sum(~isnan(tmpdat), 1);
    [~, varidx]     = max(nrdatapoints);
    varname         = dat.Properties.VariableNames{varnames(varidx)};
    
    % scatterplot
    s = scatter(dat.repetition, dat.(varname), 50, dat.session, '.');
    box off;
    xlabel('p(repeat)'); ylabel('v ~ prevresp');
    title(datasetnames{d});
    lsline;
    
    % layout
    % axis tight;
    axis square;
    xlim([0.4 0.6]); ylim([-0.3 0.3]);
    set(gca, 'ytick', [-0.3:0.1:0.3]);
    
end

print(gcf, '-dpdf', '~/Data/serialHDDM/fig2_driftrate.pdf');

