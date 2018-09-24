function kostis_plotOU_BIC

% ============================================ %
% DIC COMPARISON BETWEEN RAMPING MODELS
% ============================================ %
global mypath 
results = readtable('/Users/urai/Data/HDDM/summary/Anke_MEG_transition/allindividualresults_kostis.csv');
mdldic = mean([results.ouK_sp_bic results.ouK_input_bic results.ouK_lambda_bic results.ouK_vanilla_bic]);

close all;
subplot(4,5,1);
axis square; hold on;

cols2 = cbrewer('qual', 'Dark2', 8);
colors = cols2([5 3 4], :);

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
        text(i, mdldic(i) - 0.02*range(get(gca, 'ylim')), ...
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

ylabel({'\DeltaBIC from O-U'; 'without history'}, 'interpreter', 'tex');
set(gca, 'xtick', 1:3, 'xticklabel', {'offset bias', 'input bias', 'leak bias'}, 'xticklabelrotation', -30);

drawnow; tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/BIC_OU.pdf'));

