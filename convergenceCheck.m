function convergenceCheck(dat, folder)

% ============================================ %
% SHOW CONVERGENCE CHECKS 
% on all parameters of a HDDM table
% ============================================ %


for v = 2:10, %length(dat.Properties.VariableNames),
    
    thisdat = dat{:, v};
    
    figure;
    subplot(331); 
    plot(thisdat); 
    axis tight; box off; offsetAxes;
    xlabel('Samples'); ylabel('Parameter estimate');
    
    subplot(332); 
    [r,lags] = xcorr(thisdat, 1000, 'unbiased');
    plot(lags, r);
    axis tight; box off; offsetAxes;
    xlabel('Lags'); ylabel('Autocorrelation');
    
    subplot(333);
    histogram(thisdat, 'edgecolor', 'none');
    axis tight; box off; offsetAxes;
    xlabel('Parameter estimate'); ylabel('Frequency');
    
    suplabel(dat.Properties.VariableNames{v}, 'x');
   % export_fig(gcf, sprintf('%s/var%d.pdf', folder), '-append');

end
