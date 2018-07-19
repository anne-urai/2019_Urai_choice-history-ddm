function multiplicative_vbias

% Code to fit the history-dependent drift diffusion models described in
% Urai AE, Gee JW de, Donner TH (2018) Choice history biases subsequent evidence accumulation. bioRxiv:251595
%
% MIT License
% Copyright (c) Anne Urai, 2018
% anne.urai@gmail.com

close all; clc;
addpath(genpath('~/code/Tools'));

global mypath datasets datasetnames colors

% ============================================ %
% ONE LARGE PLOT WITH PANEL FOR EACH DATASET
% ============================================ %

for d = [4 5] % only NatComm and Anke_MEG_neutral, with varying coherence level
    disp(datasets{d});
    figure;
    
    
    results = readtable(sprintf('%s/summary/%s/allindividualresults.csv', mypath, datasets{d}));
    results = results(results.session == 0, :);
    
    allresults = struct(); alltitles = {};
    switch d
        case 5
            cohs = {'00625'  '0125' '025' '05' '1' '2' '3'};
            cohlevels = [0.625 1.25 2.5 5 10 20 30];
        case 4
            cohs = {'0' '03' '09' '27' '81'};
            cohlevels = [0 3 9 27 81];
            
    end
    
    scattercols  = cbrewer('seq', 'PuBuGn', numel(unique(cohs)) + 5);
    scattercols  = scattercols([3:end-4 end], :);
    colormap(scattercols);
    
    for c = 1:length(cohs),
        % dc_0_3_1__stimcodingdczprevrespmultiplicative
        allresults(c).v_prevresp        = results.(sprintf('dc_0_%s_1__stimcodingdczprevrespmultiplicative', cohs{c})) - ...
            results.(sprintf('dc_0_%s_2__stimcodingdczprevrespmultiplicative', cohs{c}));
        thiscoh = num2str(cohlevels(c));
        if d == 4,
            % thiscoh = num2str(cohlevels(c) * 100);
        end
        allresults(c).criterionshift    = results.(['repetition_c' regexprep(thiscoh, '\.', '\_')]);
        allresults(c).subjnr            = results.subjnr;
        allresults(c).marker 			= 'o';
        allresults(c).meancolor 		= 0.9*scattercols(c, :);
        allresults(c).scattercolor	 	= scattercols(c, :);
    end
    
    %% FLIP AROUND THE BIAS TERM FOR ALTERNATORS, THEN PLOT ITS MAGNITUDE AS A FUNCTION OF COHERENCE LEVELS
    
    subplot(4,4,1); hold on;
    vbias_all = [allresults(:).v_prevresp]; % nsj x ncoh
    vbias_all(results.repetition < 0.5, :) =  - vbias_all(results.repetition < 0.5, :);
    
    %% start with baseline: vbias (in direction of bias) for the single-vbias model
    vbias_single = results.dc_1__stimcodingdczprevresp - results.dc_2__stimcodingdczprevresp;
    vbias_single(results.repetition < 0.5, :) =  - vbias_single(results.repetition < 0.5, :);
    
    plot([1 length(cohlevels)], [nanmean(vbias_single) nanmean(vbias_single)], '-', 'color', [0 0 0]);
    s1 =  scatter(1:length(unique(cohlevels)), nanmean(vbias_single)*ones(1, length(unique(cohlevels))), 25, unique(cohlevels), 'o', 'filled');
    s1.MarkerEdgeColor = 'k';
    
    %% then the model with difficulty-dependent vbias
    
    plot(1:length(unique(cohlevels)),   nanmean(vbias_all), '-', 'color',[0 0 0]);
    s2 = scatter(1:length(unique(cohlevels)),   nanmean(vbias_all), 25, unique(cohlevels), 's', 'filled');
    s2.MarkerEdgeColor = 'k';
    
    % plot(1:length(cohlevels), (vbias_all), '-', 'color', [0.8 0.8 0.8], 'linewidth', 0.5);
    %boundedline(1:length(cohlevels), nanmean(vbias_all), nanstd(vbias_all) ./ sqrt(length(results.repetition)));
    set(gca, 'xtick', 1:length(cohlevels), 'xticklabel', cohlevels);
    
    box off; axis tight; axis square;
    set(gca, 'ytick', unique([0 get(gca, 'ytick')])); ylim([min([0 min(get(gca, 'ylim'))]) max(get(gca, 'ylim'))]);
    title(datasetnames{d});
    
    switch d
        case 4
            ylim([0 0.4]); set(gca, 'ytick', 0:0.1:0.4);
        case 5
            ylim([0 0.3]); set(gca, 'ytick', 0:0.1:0.3);
    end
    xlim([min([0.5 get(gca, 'xlim')]) max(get(gca, 'xlim'))]);
    offsetAxes;
    
    ylabel({'History shift in v_{bias}' 'towards preference'});
    
    %% ADD BETA
    [b_all, ~, stats] = glmfit(nanmean(vbias_all)', 1:length(cohlevels));
    disp([b_all(2) stats.p(2)]);
    
    [r, pval] = corr(nanmean(vbias_all)', transpose(1:length(cohlevels)))
    
    if pval < 0.001,
        text(length(cohlevels)*0.8, mean(get(gca, 'ylim')), sprintf('r %0.2f, ***', r), 'fontsize', 6);
        
    elseif pval < 0.01,
        text(length(cohlevels)*0.8, mean(get(gca, 'ylim')), sprintf('r %0.2f, **', r), 'fontsize', 6);
        
    elseif pval < 0.05,
        text(1.5, mean(get(gca, 'ylim'))*0.3, sprintf('r %0.2f, *', r), 'fontsize', 5);
    end
    
    switch d
        case 4
            xlabel('Coherence (%)');
            
        case 5
            xlabel('\Delta coherence (%)');
            set(gca, 'xticklabelrotation', -30);
    end
    
    
    tightfig;
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/vbias_per_coherence_d%d.pdf', d));
    %
    %     % for every subject, get a beta weight
    %     b = nan(size(vbias_all, 1), 2);
    %     for sj = 1:size(vbias_all, 1),
    %         b(sj, :) = glmfit(vbias_all(sj, :)', 1:length(cohlevels));
    %
    %         r = corr(vbias_all(sj, :)', transpose(1:length(cohlevels)), 'type', 'pearson');
    %         b(sj, 2) = r;
    %     end
    %     figure;
    %     subplot(221); plotBetasSwarm(b(results.repetition < 0.5, 2)); title('Alternators');
    %     subplot(222); plotBetasSwarm(b(results.repetition > 0.5, 2));  title('Repeaters');
    %     subplot(223); plotBetasSwarm(b(:, 2));  title('All');
    %
    %
    %     %% make a nice looking, small subplot
    %     figure;
    %     subplot(381);
    %     scatter(ones(size(b(:, 2))), b(:, 2), 10, [0.5 0.5 0.5], 'jitter', 'on', 'jitteramount', 0.1);
    %     hold on;
    %     plot([0.9 1.1], [nanmean(b(:, 2)) nanmean(b(:, 2)) ], 'k');
    %     % xlim([0.85 1.105]);
    %
    %     p = permtest(b(:, 2));
    %     mysigstar(gca, 1, max(get(gca, 'ylim')), p);
    %
    %     set(gca, 'xtick', [0.95 1.05], 'xticklabel', []);
    %     ylabel('\rho');
    %     %axis tight;
    %     offsetAxes;
    %     tightfig;
    %     print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/vbias_per_coherence_beta_ d%d.pdf', d));
    %
end

for d = [4 5],
    
    disp(datasets{d});
    
    
    results = readtable(sprintf('%s/summary/%s/allindividualresults.csv', mypath, datasets{d}));
    results = results(results.session == 0, :);
    
    
    
end

