function kostis_plotOU

global mypath datasets datasetnames

%% compare correlation coefficients
results = readtable(sprintf('%s/summary/%s/allindividualresults_ou.csv', mypath, 'Anke_MEG_Neutral'));
corrplot(results, {'repetition_alldata', 'ouK_input_inputbias', 'ouK_lambda_lambdabias', 'ouK_sp_spbias'});
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/OU_scatter.pdf'));

%% separately, test correlation coefficients
    
y = [results.ouK_input_inputbias results.ouK_lambda_lambdabias];
x = results.repetition_alldata;
close all;

for yidx = 1:2,
    
    thisax = subplot(4,4,yidx); % plot in two subplots
    hold on; axis square; box on;
    
    % correlate the two
    [rho, pval] = corr(x, y(:, yidx), 'type', 'pearson');
    if pval < 0.05,
        % create my own regression line, or see lsline
        b = regress(y(:, yidx), [ones(length(x), 1) x]);
        plot(x, b(1) + x*b(2), 'k-');
    end
    
    scatter(x, y(:, yidx), 10, 'o', 'filled'); % axis tight; box off;
    if pval < 0.001,
        title(sprintf('r = %.2f, p < 0.001', rho), 'fontweight', 'normal');
    else
        title(sprintf('r = %.2f, p = %.2f', rho, pval), 'fontweight', 'normal');
    end
    
    switch yidx
        case 1
            ylabel('Input bias');
        case 2
            ylabel('Leak bias');
    end
    
    % put the y axis on the right, and make sure the label is rotated and
    % moved into the right position
    if yidx == 2,
        ax = gca;
        ax.YLabel.Rotation = 270;
        ax.YAxisLocation = 'right';
        axpos = ax.YLabel.Position;
        axpos(1) = axpos(1) + 2;
        ax.YLabel.Position = axpos;
    end
end

% move the right subplot closer towards the left one
spos = get(gca, 'position');
spos(1) = 0.8*spos(1);
set(gca, 'position', spos);

% shared x label
s       = suplabel('P(repeat)', 'x');
spos    = s.Position;
spos(2) = spos(2) + 0.05; % move  up
s.Position = spos;

% test if those two correlations are different using Steiger's test
[rddiff,cilohi,p] = rddiffci(corr(x, y(:, 1)), corr(x, y(:, 2)), ...
    corr(y(:, 1), y(:, 2)), length(x), 0.05);

% plot on top
if p < 0.001,
    [a, h] = suplabel(sprintf('delta r = %.3f, p < 0.001', rddiff), 't');
else
    [a, h] = suplabel(sprintf('delta r = %.3f, p = %.3f', rddiff, p), 't');
end
set(h, 'fontweight', 'normal');
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/OU_correlations.pdf'));


end