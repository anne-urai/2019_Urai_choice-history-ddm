function correlate_z_vbias

% Code to fit the history-dependent drift diffusion models as described in
% Urai AE, de Gee JW, Tsetsos K, Donner TH (2019) Choice history biases subsequent evidence accumulation. eLife, in press.
%
% MIT License
% Copyright (c) Anne Urai, 2019
% anne.urai@gmail.com

close all; clc;
addpath(genpath('~/code/Tools'));

global mypath datasets datasetnames colors

% ============================================ %
% ONE LARGE PLOT WITH PANEL FOR EACH DATASET
% ============================================ %

for d = length(datasets):-1:1
    
    results = readtable(sprintf('%s/summary/%s/allindividualresults.csv', mypath, datasets{d}));
    results = results(results.session == 0, :);
    
    whichmdls = 'stimcoding';
    try
        % use the stimcoding difference
        results.z_prevresp = ...
            results.(['z_1__' whichmdls 'dczprevresp']) - results.(['z_2__' whichmdls 'dczprevresp']);
        results.v_prevresp = ...
            results.(['dc_1__' whichmdls 'dczprevresp']) - results.(['dc_2__' whichmdls 'dczprevresp']);
    catch
        results.z_prevresp = ...
            results.(['z_1_0__' whichmdls 'dczprevresp']) - results.(['z_2_0__' whichmdls 'dczprevresp']);
        results.v_prevresp = ...
            results.(['dc_1_0__' whichmdls 'dczprevresp']) - results.(['dc_2_0__' whichmdls 'dczprevresp']);
    end
    
    close all;
    subplot(441); hold on;
    
    % COMPUTE THE CORRELATION COEFFICIENT
    [rho, pval] = corr( results.v_prevresp,  results.z_prevresp, ...
        'type', 'spearman', 'rows', 'complete');
    
    if pval < 0.05,
        
        % 18 SEPTEMBER 2018 - SWITCH TO PCA CORRELATION https://elifesciences.org/articles/00638
        b = deming(results.v_prevresp, results.z_prevresp);
        xrangeextra = 0.15*range(results.v_prevresp);
        xrange = linspace(min(results.v_prevresp)- xrangeextra, ...
            max(results.v_prevresp)+xrangeextra, 100);
        yrange = polyval(fliplr(b'), xrange);
        
        % NOW PLOT
        l = plot(xrange, yrange);
        l.Color = 'k';
        l.LineWidth = 0.5;
        l.LineStyle = '-';
    end
    
    allresults.rho(d)   = rho;
    allresults.pval(d)  = pval;
    allresults.bf(d)    = corrbf(rho, numel(results.z_prevresp));
    
    % PLOT ALL DATAPOINTS IN SPECIFIC COLOR
    s  = scatter(results.v_prevresp,  results.z_prevresp, 7, [0.5 0.5 0.5], 'o');
    
    ylabel('History shift in z', 'color', colors(1, :));
    xlabel('History shift in v_{bias}', 'color', colors(2, :));
    axis tight; offsetAxes;
    
    fprintf('%s, rho = %.3f, pval = %.3f \n', datasets{d}, rho, pval);
    
    title(datasetnames{d});
    tightfig;
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/zbias_vs_vbias_d%d.pdf', d));
    
end

allresults
bf = prod([allresults(:).bf])

end
