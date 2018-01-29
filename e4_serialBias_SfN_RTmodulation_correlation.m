function e4_serialBias_SfN_RTmodulation_correlation

% Code to fit the history-dependent drift diffusion models described in
% Urai AE, Gee JW de, Donner TH (2018) Choice history biases subsequent evidence accumulation. bioRxiv:251595
%
% MIT License
% Copyright (c) Anne Urai, 2018
% anne.urai@gmail.com

addpath(genpath('~/code/Tools'));
warning off; close all; clear;
global datasets datasetnames

colors = linspecer(5); % red blue green

% ========================================== %
% MODULATION OF SERIAL CHOICE BIAS
% ========================================== %

for d = 1:length(datasets),
    
    close all;
    dat = readtable(sprintf('~/Data/HDDM/summary/%s/allindividualresults.csv', ...
        datasets{d}));
    
    subplot(4,4,1); hold on;
    plotScatter(dat.v_prevresp__regressdczprevrespstimrtpupil, ...
        dat.v_prevrespprevpupil__regressdczprevrespstimrtpupil, 0.2, colors(4,:));
    xlabel('v ~ prevresp');
    ylabel('v ~ prevresp * prevpupil');
    axis square;
    %  axis tight;
    
    % also pupil
    subplot(4,4,2); hold on;
    plotScatter(dat.v_prevresp__regressdczprevrespstimrtpupil, ...
        dat.v_prevrespprevrt__regressdczprevrespstimrtpupil, 0.2, colors(4, :));
    xlabel('v ~ prevresp');
    ylabel('v ~ prevresp * prevrt');
    axis square;
    %axis tight;
    
    subplot(4,4,5); hold on;
    plotScatter(dat.z_prevresp__regressdczprevrespstimrtpupil, ...
        dat.z_prevrespprevpupil__regressdczprevrespstimrtpupil, 0.2, colors(5, :));
    xlabel('z ~ prevresp');
    ylabel('z ~ prevresp * prevpupil');
    axis square;
    %  axis tight;
    
    % also pupil
    subplot(4,4,6); hold on;
    plotScatter(dat.z_prevresp__regressdczprevrespstimrtpupil, ...
        dat.z_prevrespprevrt__regressdczprevrespstimrtpupil, 0.2, colors(5, :));
    xlabel('z ~ prevresp');
    ylabel('z ~ prevresp * prevrt');
    axis square;
    %axis tight;
    
    
    try
        ss = suplabel(cat(2, datasetnames{d}{1}, ' - ', datasetnames{d}{2}), 't');
    catch
        ss = suplabel(datasetnames{d}{1}, 't');
    end
    set(ss, 'fontweight', 'normal');
    ss.FontWeight = 'normal';
    ss.Position(2) = ss.Position(2) - 0.007;
    
    tightfig;
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/fig3_RTmodulation_correlation_d%d.pdf', d));
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
plot(xlims, [0.5 0.5], 'color', [0.5 0.5 0.5], 'linewidth', 0.2);

scatter(x,y, 15, ...
    'LineWidth', 0.001, ...
    'markeredgecolor', 'w', 'markerfacecolor', plotColor);

txt = {sprintf('r_{%d} = %.3f', length(find(~isnan(y)))-2, rho) sprintf('p = %.3f', pval)};
if pval < 0.001,
    txt = {sprintf('r_{%d} = %.3f', length(find(~isnan(y)))-2,rho) sprintf('p < 0.001')};
end
tt = text(min(get(gca, 'xlim')) + legendWhere*(range(get(gca, 'xlim'))), ...
    min(get(gca, 'ylim')) + 0.2*(range(get(gca, 'ylim'))), ...
    txt, 'fontsize', 5);
set(gca, 'color', 'none');

end
