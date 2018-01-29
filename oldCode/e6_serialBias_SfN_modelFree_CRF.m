function e6_serialBias_SfN_modelFree_CRF
% ========================================== %
% conditional response functions from White & Poldrack
% ========================================== %

addpath(genpath('~/code/Tools'));
warning off; close all; clear;
global datasets datasetnames

% ========================================== %
% MODELFREE MEASURE OF BIAS
% RT DISTRIBUTIONS
% ========================================== %

clf;
for d = 1:length(datasets),
    
    % load data
    close all;
    clearvars -except d datasets cnt datasetnames whichSJ ylims;
    csvfile = dir(sprintf('~/Data/HDDM/%s/*.csv', datasets{d}));
    csvfile = csvfile(arrayfun(@(x) ~strcmp(x.name(1),'.'), csvfile)); % remove hidden files
    alldata = readtable(sprintf('~/Data/HDDM/%s/%s', datasets{d}, csvfile.name));
    
    % define repeaters and alternators based on dc
    dat = readtable(sprintf('~/Data/HDDM/summary/%s/allindividualresults.csv', datasets{d}));
    dat = dat(dat.session == 0, :);
    repeaters = (dat.repetition > 0.5);
    
    % recode into repeat and alternate
    alldata.repeat = zeros(size(alldata.response));
    alldata.repeat(alldata.response == (alldata.prevresp > 0)) = 1;
    
    % divide RT into quantiles for each subject
    qntls = [.2, .4, .6, .8, .95]; % White & Poldrack
    qntls = [.1, .3, .5, .7, .9, 1]; % Leite & Ratcliff
    
    discretizeRTs = @(x) {discretize(x, quantile(x, [0, qntls]))};
    rtbins = splitapply(discretizeRTs, alldata.rt, findgroups(alldata.subj_idx));
    alldata.rtbins = cat(1, rtbins{:});
    % scatter(1:length(alldata.rt), alldata.rt, 10, alldata.rtbins, 'filled');
    
    % get RT quantiles for choices that are in line with or against the bias
    [gr, sjidx, rtbins] = findgroups(alldata.subj_idx, alldata.rtbins);
    cpres        = array2table([sjidx, rtbins], 'variablenames', {'subj_idx', 'rtbin'});
    cpres.choice = splitapply(@nanmean, alldata.repeat, gr);
    
    % make into a subjects by rtbin matrix
    mat = unstack(cpres, 'choice', 'rtbin');
    mat = mat{:, 2:end}; % remove the last one, only has some weird tail
    % mat = mat - 0.5;
    
    % plot
    subplot(441); hold on;
    colors = flipud(cbrewer('div', 'PuOr', 5)); % repeaters in purple, alternators in orange
    plot(qntls, 0.5*ones(size(qntls)), 'color', [0.5 0.5 0.5]);
    
    % repeaters
    boundedline(qntls, nanmean(mat(repeaters, :), 1), nanstd(mat(repeaters, :), [], 1) ./ sqrt(sum(repeaters)), 'cmap', colors(1, :));
    % plot(qntls, nanmean(mat(repeaters, :), 1), 'color', colors(1, :));
    
    [h, pval] = ttest(mat(repeaters, :), 0.5); % test against zero
    plot(qntls((h > 0)), nanmean(mat(repeaters, (h > 0)), 1), '.', 'markersize', 10, 'color', colors(1, :));
    
    % alternators
    if any(~repeaters),
        boundedline(qntls, nanmean(mat(~repeaters, :), 1), nanstd(mat(~repeaters, :), [], 1) ./ sqrt(sum(~repeaters)), 'cmap', colors(end, :));
        % plot(qntls, nanmean(mat(~repeaters, :), 1), 'color', colors(end, :));
        
        [h, pval] = ttest(mat(~repeaters, :), 0.5); % test against zero
        hold on; plot(qntls((h > 0)), nanmean(mat(~repeaters, (h > 0)), 1), '.', 'markersize', 10, 'color', colors(end, :));
        
        % also test the two against each other
        [h, pval] = ttest2(mat(repeaters, :), mat(~repeaters, :));
        hold on; plot(qntls((h > 0)), 0.5*ones(size(qntls((h > 0)))), '.', 'markersize', 12, 'color', 'k');
        
    end
    
    axis tight; box off;
    set(gca, 'xtick', qntls);
    % set(gca, 'ytick', [0 0.1]); ylim([0 0.1]);
    
    axis square;  offsetAxes;
    xlabel('RT (quantiles)'); ylabel('P(repeat)');
    title(datasetnames{d}{1});
    
    tightfig;
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/CRF_d%d.pdf', d));
    
end

end
