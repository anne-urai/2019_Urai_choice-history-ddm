function e2_serialBias_SfN_SanityChecks

global mypath datasets datasetnames
addpath(genpath('~/code/Tools'));
warning off; close all;

cmap = viridis(256);
colormap(cmap);

% ========================================== %
% DOES DRIFT RATE CORRELATE WITH D'?
% ========================================== %

for d = 1:length(datasets),
    
    close all; subplot(4,4,1);
    dat = readtable(sprintf('%s/summary/%s/allindividualresults.csv', mypath, datasets{d}));
    
    if any(strcmp(dat.Properties.VariableNames', 'v__stimcodingnohist') > 0),
        % one coherence level, grey markers
        xlabel('d'''); ylabel('Drift rate (v)');
        alldprime = dat.dprime;
        alldrift = dat.v__stimcodingnohist;
        g = gscatter(alldprime(:), alldrift(:), 1:length(alldprime), ...
            [0.3 0.3 0.3], [], 2, [], 0);
        for gi = 1:length(g),
            g(gi).Marker          = 'o';
            g(gi).MarkerEdgeColor = g(gi).Color;
            g(gi).MarkerFaceColor = 'none';
            g(gi).LineWidth       = 0.01;
        end
    else
        
        % ========================================== %
        % Anke's data - drift rate as a function of coherence
        % ========================================== %
        
        vars    = dat.Properties.VariableNames';
        cohvars = vars(~cellfun(@isempty, strfind(vars, 'dprime_c')));
        alldprime = dat{dat.session == 0, cohvars};
        
        driftvars = regexp(vars, 'v_c\w+__stimcodingnohist$', 'match');
        driftvars = vars((~cellfun(@isempty, driftvars)));
        alldrift  = dat{dat.session == 0, driftvars};
        
        if isempty(alldrift),
            driftvars = regexp(vars, 'v_0_\w+__stimcodingnohist$', 'match');
            driftvars = vars((~cellfun(@isempty, driftvars)));
            alldrift  = dat{dat.session == 0, driftvars};
        end
        
        % which coherences were used?
        cohvars = regexprep(cohvars, '_', '.');
        cohs    = cellfun(@sscanf, cohvars, repmat({'dprime.c%f'}, length(cohvars), 1));
        allcohs = repmat(cohs', size(alldrift, 1), 1);
        colors  = cbrewer('seq', 'PuBuGn', numel(unique(allcohs(:))) + 5);
        colors  = colors([3:end-4 end], :);
        g       = gscatter(alldprime(:), alldrift(:), allcohs(:), colors, [], 2, [], 0);
        
        box off;
        for gi = 1:length(g),
            g(gi).Marker = 'o';
            g(gi).MarkerEdgeColor = g(gi).Color;
            g(gi).MarkerFaceColor = 'none';
            g(gi).LineWidth = 0.01;
        end
    end
    
    %% add correlation
    % layout
    axis tight; axis square;
    
    switch datasets{d}
        case 'JW_yesno'
           ylim([0.4 1.09]); 
    end
    xlim([0 ceil(max(get(gca, 'xlim')))]);
    xlabel('d''');
    if d == 1,
        ylabel('Drift rate (v)');
    end
    
    if any(strcmp(dat.Properties.VariableNames', 'v__stimcodingnohist') > 0),
        [rho, pval] = corr(alldprime(:), alldrift(:), 'rows', 'complete');
    else 
        [coef,pval] = partialcorr([alldprime(:), alldrift(:)], allcohs(:), 'rows', 'complete');
        rho = coef(1,2); pval = pval(1,2);
    end
    
    txt = {sprintf('r_{%d} = %.3f',  length(find(~isnan(alldrift(:))))-2, rho) sprintf('p = %.3f', pval)};
    if pval < 0.001,
        txt = {sprintf('r_{%d} = %.3f',  length(find(~isnan(alldrift(:))))-2, rho) sprintf('p < 0.001')};
    end
    tt = text(min(get(gca, 'xlim')) + 0.04*(range(get(gca, 'xlim'))), ...
        min(get(gca, 'ylim')) + 0.8*(range(get(gca, 'ylim'))), ...
        txt, 'fontsize', 5);
    set(gca, 'color', 'none');
    
    % lsline manually
    p = polyfit(alldprime(~isnan(alldrift) & ~isnan(alldprime)), ...
        alldrift(~isnan(alldrift(:)) & ~isnan(alldprime(:))), 1);
    l = refline(p(1), p(2));
    l.Color = 'k';
    l.LineWidth = 0.5;
    if pval < 0.05,
        l.LineStyle = '-';
    else
        l.LineStyle = ':';
    end
    
    offsetAxes; box off;
    title(datasetnames{d});
    
    tightfig;
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/figure1b_HDDM_driftrate_d%d.pdf',d));
    
    % ========================================== %
    % colorbar
    % ========================================== %
    
    if ~any(strcmp(dat.Properties.VariableNames', 'v__stimcodingnohist') > 0),
        
        close all;
        subplot(441);
        colormap(colors);
        imagesc(log(allcohs+1));
        handles = colorbar('southoutside');
        set(gca, 'Xscale', 'log');
        
        % ==================================================================
        % make colorbar look pretty
        % ==================================================================
        drawnow;
        handles.TickDirection = 'out';
        handles.Box = 'off';
        
        % find the right place for tickmarks
        if d == 2,
            handles.Ticks = linspace(0.75, 3.25, 7);
            handles.TickLabels = {'0.6' '1.2' '2.5' '5' '10' '20' '30'};
        elseif d == 4,
            handles.Ticks = linspace(0.5, 3.75, 6);
            handles.TickLabels = {'0' '5' '10' '20' '40', '60'};
        elseif d == 9,  
            handles.Ticks = linspace(0.35, 4.2, 10);
            handles.TickLabels = {'0' '3' '5' '9' '10', '20', '27', '40', '60', '81'};
        end
        drawnow;
        
        % get original axes
        hAllAxes = findobj(gcf,'type','axes'); axpos = {};
        axpos = hAllAxes.Position;
        
        % make colorbar thinner
        cpos = handles.Position;
        cpos(3) = 0.75*cpos(3);
        handles.Position = cpos;
        
        drawnow;
        % restore axis pos
        set(hAllAxes, 'Position', axpos);
        drawnow;
        handles.Label.String = '% coherence';
        handles.FontSize = 5;
        
        % export_fig(handles, sprintf('~/Data/serialHDDM/figure1b_legend_d%d.pdf',d));
        axis off;
        % tightfig;
        print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/figure1b_legend_d%d.pdf',d));
        
    end
    
end
