
% ============================================ %
% DIC COMPARISON BETWEEN RAMPING MODELS
% ============================================ %
global mypath colors
colors(3, :) = mean(colors([1 2], :));

results = readtable('/Users/urai/Data/HDDM/summary/Anke_MEG_transition/allindividualresults_kostis.csv');
mdldic = mean([results.ddmK_z_bic results.ddmK_dc_bic  results.ddmK_dcz_bic results.ddmK_vanilla_bic]);

close all;
subplot(4,5,1);
axis square; hold on;

% cols1 = cbrewer('qual', 'Set1', 8);
% cols2 = cbrewer('qual', 'Dark2', 8);
% colors = [cols1(2, :); cols2(6, :); nanmean([cols1(2, :); cols2(6, :)])];

%everything relative to the full model
mdldic = bsxfun(@minus, mdldic, mdldic(end));
mdldic = mdldic(1:end-1);
[~, bestMdl] = min(mdldic);

for i = 1:length(mdldic),
    b = bar(i, mdldic(i), 'facecolor', colors(i, :), 'barwidth', 0.6, 'BaseValue', 0, ...
        'edgecolor', 'none');
end

%# Add a text string above/below each bin
for i = 1:length(mdldic),
    if mdldic(i) < 0,
        text(i, mdldic(i) + 0.12*range(get(gca, 'ylim')), ...
            num2str(round(mdldic(i))), ...
            'VerticalAlignment', 'top', 'FontSize', 4, 'horizontalalignment', 'center', 'color', 'w');
    elseif mdldic(i) > 0,
        text(i, mdldic(i) - 0.12*range(get(gca, 'ylim')), ...
            num2str(round(mdldic(i))), ...
            'VerticalAlignment', 'top', 'FontSize', 4, 'horizontalalignment', 'center', 'color', 'w');
    end
end
axis square; axis tight;
xlim([0.5 length(mdldic)+0.5]);
offsetAxes; box off;
set(gca, 'color', 'none');
set(gca, 'xcolor', 'k', 'ycolor', 'k');
%title(datasetnames{d});

ylabel({'\DeltaBIC from DDM'; 'without history'}, 'interpreter', 'tex');
set(gca, 'xtick', 1:3, 'xticklabel', {'z', 'v_{bias}', 'z+v_{bias}'}, 'xticklabelrotation', -30);

drawnow; tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/BIC_DDM.pdf'));

