function correlations_lags_regression

global mypath datasets datasetnames
addpath(genpath('~/code/Tools'));
warning off; close all;
colors = cbrewer('qual', 'Set2', length(datasets));

set(groot, 'defaultaxesfontsize', 6, 'defaultaxestitlefontsizemultiplier', 1.1, ...
    'defaultaxeslabelfontsizemultiplier', 1, ...
    'defaultaxestitlefontweight', 'bold', ...
    'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters', ...
    'DefaultAxesBox', 'off', ...
    'DefaultAxesTickLength', [0.02 0.05], 'defaultaxestickdir', 'out', 'DefaultAxesTickDirMode', 'manual', ...
    'defaultfigurecolormap', [1 1 1], 'defaultTextInterpreter','tex', ...
    'DefaultFigureWindowStyle','normal');

% ==================================== %
% CORRELATIONS PER LAG
% ==================================== %

close all;
for d = 1:length(datasets),
    dat = readtable(sprintf('%s/summary/%s/allindividualresults.csv', mypath, datasets{d}));
    
    % compute correlation with p_repeat
    for l = 1:3,
        if l == 1,
            [mat_z.r(d, l), mat_z.ci(d,l,:), mat_z.pval(d,l)] = ...
                spearmans(dat.z_prevresp__regressdczlag3, dat.repetition);
            [mat_dc.r(d, l), mat_dc.ci(d,l,:), mat_dc.pval(d,l)] = ...
                spearmans(dat.v_prevresp__regressdczlag3, dat.repetition);
        else
            [mat_z.r(d, l), mat_z.ci(d,l,:), mat_z.pval(d,l)] = ...
                spearmans(dat.(['z_prev' num2str(l) 'resp__regressdczlag3']), ...
                dat.(['repetition' num2str(l)]));
            [mat_dc.r(d, l), mat_dc.ci(d,l,:), mat_dc.pval(d,l)] = ...
                spearmans(dat.(['v_prev' num2str(l) 'resp__regressdczlag3']), ...
                dat.(['repetition' num2str(l)]));
        end
    end
end

% Z CORRELATION KERNELS
close all;
subplot(331); hold on;
plot([1 3], [0 0], 'k', 'linewidth', 0.5);
for d = 1:length(datasets),
    plot(1:3, mat_z.r(d, :), 'color', colors(d, :), 'linewidth', 1);
    h = (mat_z.r(d,:) < 0.05);
    if any(h>0),
        plot(find(h==1), mat_z.r(d, (h==1)), '.', 'markeredgecolor', colors(d, :), ...
        'markerfacecolor', colors(d,:), 'markersize', 5);
    end
end

plot(1:3, nanmean(mat_z.r), 'k', 'linewidth', 1);
[h, pval] = ttest(mat_z.r);
if any(h>0),
    plot(find(h==1), nanmean(mat_z.r(:, (h==1))), 'ok', 'markeredgecolor', 'w', 'markerfacecolor', 'k', 'markersize', 4);
end

ylabel({'Correlation \rho, P(repeat)' 'with z ~ prevresp'})
ylim([-0.5 1]);
set(gca, 'xtick', 1:3, 'xcolor', 'k', 'ycolor', 'k');
xlabel('Lag (past trials)');
axis square; offsetAxes; tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/correlationkernels_z.pdf'));

% DC CORRELATION KERNELS
close all;
subplot(331); hold on;
plot([1 3], [0 0], 'k', 'linewidth', 0.5);
for d = 1:length(datasets),
    plot(1:3, mat_dc.r(d, :), 'color', colors(d, :), 'linewidth', 1);
    h = (mat_dc.pval(d,:) < 0.05);
    if any(h>0),
        plot(find(h==1), mat_dc.r(d, (h==1)), '.', 'markeredgecolor', colors(d, :), ...
        'markerfacecolor', colors(d,:), 'markersize', 5);
    end
end

plot(1:3, nanmean(mat_dc.r), 'k', 'linewidth', 1);
[h, pval] = ttest(mat_dc.r);
if any(h>0),
    plot(find(h==1), nanmean(mat_dc.r(:, (h==1))), 'ok', 'markeredgecolor', 'w', 'markerfacecolor', 'k', 'markersize', 4);
end

ylabel({'Correlation \rho, P(repeat)' 'with dc ~ prevresp'})
ylim([-0.5 1]);
set(gca, 'xtick', 1:3, 'xcolor', 'k', 'ycolor', 'k');
xlabel('Lag (past trials)');
axis square; offsetAxes; tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/correlationkernels_dc.pdf'));

end