function multiplicative_vbias

% Code to fit the history-dependent drift diffusion models described in
% Urai AE, Gee JW de, Donner TH (2018) Choice history biases subsequent evidence accumulation. bioRxiv:251595
%
% MIT License
% Copyright (c) Anne Urai, 2018
% anne.urai@gmail.com

close all; clc;
addpath(genpath('~/code/Tools'));

global mypath datasets datasetnames

for d = [4 5],
    
    disp(datasets{d});
    
        results = readtable(sprintf('%s/summary/%s/allindividualresults.csv', mypath, datasets{d}));
    results = results(results.session == 0, :);
      switch d
        case 5
            cohs = {'c0_0625'  'c1_25' 'c2_5' 'c5' 'c10' 'c20' 'c30'};
            cohlevels = [0.625 1.25 2.5 5 10 20 30];
        case 4
            cohs = {'0' '3' '9' '27' '81'};
            cohlevels = [0 3 9 27 81];
            
    end
    
    scattercols  = cbrewer('seq', 'PuBuGn', numel(unique(cohs)) + 5);
    scattercols  = scattercols([3:end-4 end], :);
    
    clear multiplicative single coherencecolors;
    for c = 1:length(cohs),
        % dc_0_3_1__stimcodingdczprevrespmultiplicative
        multiplicative(:, c) = results.(sprintf('v_%s__stimcodingdcprevrespmultiplicative', cohs{c}));
        single(:, c) = results.(sprintf('v_%s__stimcodingdczprevresp', cohs{c}));
        coherencecolors(:, c) = c*ones(size(results.(sprintf('v_%s__stimcodingdcprevresp', cohs{c}))));
    end
    
    close all; subplot(441);
    scatter(single(:), multiplicative(:), 10, coherencecolors(:));
    xlabel('Single v_{bias}');
    ylabel('Coherence-dependent v_{bias}');
    title('Drift rate (v)');
    axis tight;  axisEqual;
    colormap(scattercols);
    axis square; offsetAxes;
    
      tightfig;
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/driftRate_comparison_d%d.pdf', d));
  
end


%% ALSO SHOW BEHAVIOR: P(BIAS) AND RT
for d = [4 5],
    
    disp(datasets{d});
    
    switch d
        case 5
            cohs = {'00625'  '0125' '025' '05' '1' '2' '3'};
            cohlevels = [0.625 1.25 2.5 5 10 20 30];
        case 4
            cohs = {'0' '03' '09' '27' '81'};
            cohlevels = [0 3 9 27 81];
    end
    
    % find the original datafile
    filename = dir(sprintf('%s/%s/*.csv', mypath, datasets{d}));
    alldata  = readtable(sprintf('%s/%s/%s', mypath, datasets{d}, filename.name));
    
    % recode into repeat and alternate for the model
    alldata.repeat = zeros(size(alldata.response));
    alldata.repeat(alldata.response == (alldata.prevresp > 0)) = 1;
    
    % for each observers, compute their bias
    [gr, sjs] = findgroups(alldata.subj_idx);
    sjrep = splitapply(@nanmean, alldata.repeat, gr);
    sjrep = sjs(sjrep < 0.5); % alternators
    
    % recode into biased and unbiased choices
    alldata.biased = alldata.repeat;
    altIdx = ismember(alldata.subj_idx, sjrep);
    alldata.biased(altIdx) = double(~(alldata.biased(altIdx))); % flip
    % biased [0,1] now codes for repeat [0,1] for repeaters and switch
    % [0,1] for alternaters
    
    if numel(unique(alldata.stimulus)) > 2,
        alldata.stimulus = sign(alldata.stimulus);
    elseif isequal(unique(alldata.stimulus), [0 1]'),
        alldata.stimulus = sign(alldata.stimulus - 0.1);
    end
    alldata.coherence       = [alldata.coherence .* alldata.stimulus];
    
    % from fraction to percent
    if all(abs(alldata.coherence) < 1),
        alldata.coherence = alldata.coherence * 100;
    end
    
    clear tab;
    [gr, sj, coh]   = findgroups(alldata.subj_idx, abs(alldata.coherence));
    tab.subj_idx    = splitapply(@mean, alldata.subj_idx, gr);
    tab.coherence   = splitapply(@mean, abs(alldata.coherence), gr);
    tab = struct2table(tab);
    tab.repetition = splitapply(@mean, alldata.biased, gr);
    tab.rt = splitapply(@mean, alldata.rt, gr);
    
    sjrep = splitapply(@nanmean, alldata.repeat, findgroups(alldata.subj_idx));
    
    % plot RT
    tab2 = unstack(tab(:, [1 2 4]), 'rt', 'subj_idx');
    tab2 = tab2{:, 2:end}';
    tab2_repeat = tab2(sjrep > 0.5, :);
    tab2_alternate = tab2(sjrep < 0.5, :);
    
    close all; subplot(4,4,1); hold on;
    errorbar(1:length(cohlevels), nanmean(tab2), ...
        nanstd(tab2) ./ sqrt(sum(~isnan(tab2))),'ok-', 'capsize', 0, 'markerfacecolor', 'w');
    %errorbar(1:length(cohlevels), nanmean(tab2_alternate), ...
    %    nanstd(tab2_alternate) ./ sqrt(sum(~isnan(tab2_alternate))),'ob-', 'capsize', 0, 'markerfacecolor', 'w');
    
    set(gca, 'xtick', 1:length(cohlevels), 'xticklabel', cohlevels);
    
    box off; axis tight; axis square;
    title(datasetnames{d})
    xlim([min([0.5 get(gca, 'xlim')]) max(get(gca, 'xlim'))]);
    offsetAxes;
    ylabel('Response time (s)');
    tightfig;
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/rt_per_coherence_d%d.pdf', d));
    
    % BIASED BEHAVIOR
    tab2 = unstack(tab(:, [1 2 3]), 'repetition', 'subj_idx');
    tab2 = tab2{:, 2:end}';
    tab2_repeat = tab2(sjrep > 0.5, :);
    tab2_alternate = tab2(sjrep < 0.5, :);
    
    close all; subplot(4,4,1); hold on;
    errorbar(1:length(cohlevels), nanmean(tab2), ...
        nanstd(tab2) ./ sqrt(sum(~isnan(tab2))),'ok-', 'capsize', 0, 'markerfacecolor', 'w');
    %errorbar(1:length(cohlevels), nanmean(tab2_alternate), ...
    %    nanstd(tab2_alternate) ./ sqrt(sum(~isnan(tab2_alternate))),'ob-', 'capsize', 0, 'markerfacecolor', 'w');
    
    set(gca, 'xtick', 1:length(cohlevels), 'xticklabel', cohlevels);
    
    box off; axis tight; axis square;
    title(datasetnames{d})
    xlim([min([0.5 get(gca, 'xlim')]) max(get(gca, 'xlim'))]);
    offsetAxes;
    ylabel('P(bias)');
    xlabel('Coherence (%)');
    tightfig;
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/repetition_per_coherence_d%d.pdf', d));
    
end


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
        allresults(c).v_prevresp        = results.(sprintf('dc_0_%s_1__stimcodingdcprevrespmultiplicative', cohs{c})) - ...
            results.(sprintf('dc_0_%s_2__stimcodingdcprevrespmultiplicative', cohs{c}));
        thiscoh = num2str(cohlevels(c));
        if d == 4,
            thiscoh = num2str(cohlevels(c) * 100);
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
    errorbar(1:length(cohlevels), nanmean(vbias_single)*ones(size(cohlevels)), ...
        nanstd(vbias_single) ./ sqrt(length(vbias_single)) *ones(size(cohlevels)),'ok-', 'capsize', 0, 'markerfacecolor', 'w');
    
    %% then the model with difficulty-dependent vbias
    errorbar(1:length(cohlevels), nanmean(vbias_all), ...
        nanstd(vbias_all) ./ sqrt(length(vbias_single)) ,'sr-', 'capsize', 0, 'markerfacecolor', 'w');
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
    
end


end

