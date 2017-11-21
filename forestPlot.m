function forestPlot(alldat)

ds = [17 10 16 15 1 2 3 7 8 9]; % from the e1b collection script
ds = length(alldat):-1:1;
close all;

global colors;
axiscolors = colors; 

% MAKE AN OVERVIEW PLOT
subplot(331); hold on;
% make a vertical line at zero
plot([0 0], [0.5 length(ds)+0.5], 'color', [0 0 0], 'linewidth', 0.5);

for d = 1:length(ds),
    
    markercolors = cbrewer('qual', 'Paired', 10);
    transitioncolors = [[0.5 0.5 0.5]; markercolors([7 9], :)];
    meancolors = [0 0 0; markercolors([8 10], :)];
    markers = {'o', 'v', '^'}; %also indicate with different markers
    
    % determine color and marker
    if ~isempty(strfind(alldat(ds(d)).datasetnames, 'Repetitive'))
        col = transitioncolors(3, :);
        mrk = markers{3};
        meancol = meancolors(3, :);
    elseif ~isempty(strfind(alldat(ds(d)).datasetnames, 'Alternating'))
        col = transitioncolors(2, :);
        mrk = markers{2};
        meancol = meancolors(2, :);
    else
        col = transitioncolors(1, :);
        mrk = markers{1};
        meancol = meancolors(1, :);
    end
    
    % start at the top
    h = ploterr(alldat(ds(d)).corrz, ...
        length(ds)-d+1,  ...
        {alldat(ds(d)).corrz_ci(1) alldat(ds(d)).corrz_ci(2)}, [], 'o', 'abshhxy', 0.2);
    set(h(1), 'marker', mrk, 'color', col, 'markerfacecolor', meancol, 'markeredgecolor', 'w', ...
        'markersize', 5 , 'linewidth', 0.5);
    set(h(2), 'color', col, 'linewidth', 1);
end

set(gca, 'ytick', 1:length(ds), 'yticklabel', fliplr({alldat(ds).datasetnames}));
xlabel('z_{bias}');
xlim([-1 1]); offsetAxes;
set(gca, 'xcolor', axiscolors(2, :));

plot(nanmean([alldat(ds).corrz]), 0.1, 'd', 'color', 'k', 'markersize', 4);
[h, pval, ci, stats] = ttest([alldat(ds).corrz]);
bf = prod([alldat(ds).bfz]);
if bf < 100,
    title(sprintf('BF_{10} < 1/100'));
elseif bf > 100,
    title(sprintf('BF_{10} > 100'));
elseif bf < 1,
    title(sprintf('BF_{10} = 1/%.2f', 1/bf));
elseif bf > 1,
    title(sprintf('BF_{10} = %.2f', bf));
end

%% NOW FOR DRIFT CRITERION
% MAKE AN OVERVIEW PLOT
sp2 = subplot(332); hold on;
% make a vertical line at zero
plot([0 0], [0.5 length(ds)+0.5], 'color', [0 0 0], 'linewidth', 0.5);

for d = 1:length(ds),
    
    markercolors = cbrewer('qual', 'Paired', 10);
    transitioncolors = [[0.5 0.5 0.5]; markercolors([7 9], :)];
    meancolors = [0 0 0; markercolors([8 10], :)];
    markers = {'o', 'v', '^'}; %also indicate with different markers
    
    % determine color and marker
    if ~isempty(strfind(alldat(ds(d)).datasetnames, 'Repetitive'))
        col = transitioncolors(3, :);
        mrk = markers{3};
        meancol = meancolors(3, :);
    elseif ~isempty(strfind(alldat(ds(d)).datasetnames, 'Alternating'))
        col = transitioncolors(2, :);
        mrk = markers{2};
        meancol = meancolors(2, :);
    else
        col = transitioncolors(1, :);
        mrk = markers{1};
        meancol = meancolors(1, :);
    end
    
    % start at the top
    h = ploterr(alldat(ds(d)).corrv, ...
        length(ds)-d+1,  ...
        {alldat(ds(d)).corrv_ci(1) alldat(ds(d)).corrv_ci(2)} , [], 'o', 'abshhxy', 0.2);
    set(h(1), 'marker', mrk, 'color', col, 'markerfacecolor', meancol, 'markeredgecolor', 'w', ...
        'markersize', 5, 'linewidth', 0.5);
    set(h(2), 'color', col, 'linewidth', 1);
end

set(gca, 'ytick', 1:length(ds), 'yticklabel', fliplr({alldat(ds).datasetnames}), 'YAxisLocation', 'right');
xlabel('v_{bias}');
set(gca, 'xcolor', axiscolors(1, :));
xlim([-1 1]); offsetAxes;

% ADD THE AVERAGE??
plot(nanmean([alldat(ds).corrv]), 0.1, 'd', 'color', 'k', 'markersize', 4);
[h, pval, ci, stats] = ttest(fisherz([alldat(ds).corrv]));

bf = prod([alldat(ds).bfv]);
if bf < 100,
    title(sprintf('BF_{10} < 1/100'));
elseif bf > 100,
    title(sprintf('BF_{10} > 100'));
elseif bf < 1,
    title(sprintf('BF_{10} = 1/%.2f', 1/bf));
elseif bf > 1,
    title(sprintf('BF_{10} = %.2f', bf));
end

% move closer together
sp2.Position(1) = sp2.Position(1) - 0.05;

%% ADD TEXT
for d = 1:length(ds),
    
    if alldat(ds(d)).pdiff < 0.001,
        txt = sprintf('\\Deltar = %.3f, ***', alldat(ds(d)).corrdiff);
    elseif alldat(ds(d)).pdiff < 0.01,
        txt = sprintf('\\Deltar = %.3f, **', alldat(ds(d)).corrdiff);
    elseif alldat(ds(d)).pdiff < 0.05,
        txt = sprintf('\\Deltar = %.3f, *', alldat(ds(d)).corrdiff);
    else
        txt = sprintf('\\Deltar = %.3f, n.s.', alldat(ds(d)).corrdiff);
    end
    text(-1.3, length(ds)-d+1, txt, ...
        'fontsize', 5);
end

% DO STATS ACROSS DATASETS!
[h, pval, ci, stats] = ttest(fisherz([alldat(ds).corrv]), fisherz([alldat(ds).corrz]));
%suplabel(sprintf('t(%d) = %.3f, p = %.4f', stats.df, stats.tstat, pval), 't');

tightfig;
% print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/forestplot.pdf'));


end

