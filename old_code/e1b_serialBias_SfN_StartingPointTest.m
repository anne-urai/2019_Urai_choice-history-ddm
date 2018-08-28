function e1b_serialBias_SfN_StartingPointTest

% Code to fit the history-dependent drift diffusion models described in
% Urai AE, Gee JW de, Donner TH (2018) Choice history biases subsequent evidence accumulation. bioRxiv:251595
%
% MIT License
% Copyright (c) Anne Urai, 2018
% anne.urai@gmail.com

clear; close all; clc;
addpath(genpath('~/code/Tools'));

global mypath datasets datasetnames

howtoplot = 'separate'
% ============================================ %
% ONE LARGE PLOT WITH PANEL FOR EACH DATASET
% ============================================ %

close all;
for d = 1:length(datasets),
    disp(datasets{d});
    
    colors = linspecer(5); % red blue green
    
    results = readtable(sprintf('%s/summary/%s/allindividualresults.csv', mypath, datasets{d}));
    results = results(results.session == 0, :);
    disp(datasets{d}); disp(numel(unique(results.subjnr)));
    
    switch howtoplot
        case 'together'
            x = results.bias;
            y1 = results.z__stimcodingnohist;
            y2 = results.dc__stimcodingnohist;
        case 'separate'
            x = results.bias;
            y1 = results.z__stimcodingnohistonlyz;
            y2 = results.dc__stimcodingnohistonlydc;
    end
    
    try
        
        close all;
        subplot(4,4,1); hold on;
        [rho1, tt1] = plotScatter(y1,x, 0.05, colors(5, :));
        xlabel('Starting point', 'interpreter', 'tex');
        ylabel('Choice bias');
        offsetAxes;
        
        sp2 = subplot(4,4,2); hold on;
        [rho2, tt2] = plotScatter(y2, x, 0.05, colors(4, :));
        xlabel('Drift criterion', 'interpreter', 'tex', 'fontweight', 'normal');
        set(gca, 'yticklabel', []);
        
        try
            % compute the difference in correlation
            rho3 = corr(results.v_prevresp__regressdczprevresp, results.z_prevresp__regressdczprevresp, ...
                'rows', 'complete', 'type', 'pearson');
            [rhodiff, ~, pval] = rddiffci(rho1,rho2,rho3,numel(~isnan(results.repetition)), 0.05);
            offsetAxes; drawnow;
            
            % move together
            sp2.Position(1) = sp2.Position(1) - 0.08;
            try
                ss = suplabel(cat(2, datasetnames{d}{1}, ' - ', datasetnames{d}{2}), 't');
            catch
                ss = suplabel(datasetnames{d}{1}, 't');
            end
            set(ss, 'fontweight', 'normal');
            ss.FontWeight = 'normal';
            ss.Position(2) = ss.Position(2) - 0.007;
            
            %% add line between the two correlation coefficients
            txt = {sprintf('\\Deltar_{%d} = %.3f, p = %.3f', length(find(~isnan(results.v_prevresp__regressdczprevresp)))-3, rhodiff, pval)};
            if pval < 0.001,
                txt = {sprintf('\\Deltar_{%d} = %.3f, p < 0.001', length(find(~isnan(results.v_prevresp__regressdczprevresp)))-3,  rhodiff)};
            end
            title(txt, 'fontweight', 'bold', 'fontsize', 6, 'horizontalalignment', 'left');
        end
        tightfig;
        print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/startingPointTest_d%d_%s.pdf',d, howtoplot));
    end
end

end

function [rho, tt] = plotScatter(x,y, legendWhere, plotColor);

plot(x,y, '.');
axis square;
axisNotSoTight;
[rho, pval] = corr(x,y, 'type', 'pearson', 'rows', 'complete');
l = lsline;

l.Color = 'k';
l.LineWidth = 0.5;
if pval < 0.05,
    l.LineStyle = '-';
else
    l.LineStyle = ':';
end

% show lines to indicate origin
xlims = [min(get(gca, 'xlim')) max(get(gca, 'xlim'))];
ylims = [min(get(gca, 'ylim')) max(get(gca, 'ylim'))];
plot([0 0], ylims, 'color', [0.5 0.5 0.5], 'linewidth', 0.2);
plot(xlims, [0 0 ], 'color', [0.5 0.5 0.5], 'linewidth', 0.2);

scatter(x,y, 15, ...
    'LineWidth', 0.001, ...
    'markeredgecolor', 'w', 'markerfacecolor', plotColor);

txt = {sprintf('r_{%d} = %.3f', length(find(~isnan(y)))-2, rho) sprintf('p = %.3f', pval)};
if pval < 0.001,
    txt = {sprintf('r_{%d} = %.3f', length(find(~isnan(y)))-2,rho) sprintf('p < 0.001')};
end
tt = text(min(get(gca, 'xlim')) + legendWhere*(range(get(gca, 'xlim'))), ...
    min(get(gca, 'ylim')) + 0.8*(range(get(gca, 'ylim'))), ...
    txt, 'fontsize', 5);
set(gca, 'color', 'none');

% also add the group mean
p = ploterr(nanmean(x), nanmean(y), nanstd(x) ./ sqrt(length(x)), ...
    nanstd(y) ./ sqrt(length(y)), '.k', 'abshhxy', 0);
set(p(1), 'markersize', 0.1); % tiny marker

end
