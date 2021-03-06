function compare_correlations_correct_error(alldat)

% Code to fit the history-dependent drift diffusion models as described in
% Urai AE, de Gee JW, Tsetsos K, Donner TH (2019) Choice history biases subsequent evidence accumulation. eLife, in press.
%
% MIT License
% Copyright (c) Anne Urai, 2019
% anne.urai@gmail.com

global datasets

vars = {'corrv', 'corrz'};
for v = 1:length(vars),
    
    % grab the two correlation coefficients from each dataset
    rho_correct = [alldat(1:2:end).(vars{v})];
    rho_error 	= [alldat(2:2:end).(vars{v})];
    n 			= [alldat(1:2:end).nsubj];
    assert(isequal(n, [alldat(2:2:end).nsubj]), 'nsubj does not match between correct and error');
    
    close all; subplot(341); hold on;
    
    % make a vertical line at zero
    plot([0 0], [0.5 length(datasets)+0.5], 'color', [0 0 0], 'linewidth', 0.5);
    cols = cbrewer('qual', 'Paired', 10);
    
    for d = 1:length(datasets),
        
        [ridiff(d),cilohi,p] = ridiffci(rho_error(d), rho_correct(d),  n(d), n(d), 0.05);
        disp(p)
        
        
        % start at the top
        h = ploterr(ridiff(d), ...
            d,  ...
            {cilohi(1) cilohi(2)}, [], 'o', 'abshhxy', 0.2);
        set(h(1), 'marker', '^', 'color', cols(9, :), 'markerfacecolor', cols(10, :), 'markeredgecolor', 'w', ...
            'markersize', 5 , 'linewidth', 0.5);
        set(h(2), 'color', cols(9, :), 'linewidth', 1);
        
        % show the p-value
        %% ADD TEXT
        if p < 0.0001,
            txt = sprintf('p < 0.0001', p);
        else
            txt = sprintf('p = %.4f', p);
        end
        text(0.3, d-0.2, txt, 'fontsize', 6);
        
    end
    
    names = {alldat(1:2:end).datasetnames};
    for n = 1:length(names),
        if numel(names{n}) > 1,
            names{n} = cat(2, names{n}{1}, ' ', names{n}{2});
        else
            names{n} = names{n}{1};
        end
    end
    
    set(gca, 'ytick', 1:length(datasets), 'yticklabel', names);
    
    % set(gca, 'ytick', 1:length(datasets), 'yticklabel', []);
    switch vars{v}
        case 'corrz'
            xlabel('\Delta\rho correct vs. error, z');
        case 'corrv'
            xlabel('\Delta\rho correct vs. error, v_{bias}');
            set(gca, 'yticklabel', []);
    end
    
    xlim([-1 1]);
    offsetAxes;
    set(gca, 'ycolor', 'k');
    
    plot(nanmean(ridiff), 0.1, 'd', 'color', 'k', 'markersize', 4);
    
    % do bayesian t-test on the difference
    [h, pval, ci, stats] = ttest(ridiff)
    bf10 = t1smpbf(stats.tstat, length(ridiff))
    t = title(sprintf('BF_{10} = %.2f', bf10), 'fontweight', 'normal', 'fontangle', 'italic');
    t.Position(2) = t.Position(2) - 1.2; % move down
    
    tightfig;
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/forestplot_correctError_%s.pdf', vars{v}));
end

end

