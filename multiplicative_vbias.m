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
   %  close all;
    
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
    
    for c = 1:length(cohs),
        % dc_0_3_1__stimcodingdczprevrespmultiplicative
        allresults(c).v_prevresp        = results.(sprintf('dc_0_%s_1__stimcodingdczprevrespmultiplicative', cohs{c})) - ...
            results.(sprintf('dc_0_%s_2__stimcodingdczprevrespmultiplicative', cohs{c}));
        thiscoh = num2str(cohlevels(c));
        allresults(c).criterionshift    = results.(['repetition_c' regexprep(thiscoh, '\.', '\_')]);
        allresults(c).subjnr            = results.subjnr;
        allresults(c).marker 			= 'o';
        allresults(c).meancolor 		= 0.9*scattercols(c, :);
        allresults(c).scattercolor	 	= scattercols(c, :);
    end
    
    
%     sp2 = subplot(3,3,1); hold on;
%     [rho2, tt2, handles] = plotScatter(allresults, 'v_prevresp', 0.05, 0);
%     set(gca, 'yticklabel', []);
%     ylabel('P(repeat)');
%     
%     % move together
%     ss = suplabel(cat(2, datasetnames{d}{1}, ' ', datasetnames{d}{2}), 't');
%     
%     set(ss, 'fontweight', 'normal');
%     ss.FontWeight = 'normal';
%     ss.Position(2) = ss.Position(2) - 0.03;
%     
%     % add colored axes after suplabel (which makes them black)
%     xlabel(sp2, 'History shift in v_{bias}');
%     set(sp2, 'xcolor', colors(2, :), 'ycolor', 'k');
%     
%     %% make another subplot with the correlation coefficients!
%     subplot(4,4,2); clear rho pval; hold on;
%     for c = 1:length(cohs),
%         [rho(c), pval(c)] = corr(allresults(c).v_prevresp, allresults(c).criterionshift, ...
%             'rows', 'complete');
%     end
%     
%     scatter(cohlevels(pval < 0.05), rho(pval < 0.05), 20, scattercols((pval < 0.05), :), 'filled');
%     scatter(cohlevels(pval > 0.05), rho(pval > 0.05), 18, scattercols((pval > 0.05), :));
%     
%     set(gca, 'xtick', cohlevels);
%     ylabel('Correlation coefficient');
%     xlabel('Coherence level');
%     offsetAxes;
% 
%     
    %% FLIP AROUND THE BIAS TERM FOR ALTERNATORS, THEN PLOT ITS MAGNITUDE AS A FUNCTION OF COHERENCE LEVELS
    
    subplot(3,3,d);
    vbias_all = [allresults(:).v_prevresp]; % nsj x ncoh
    vbias_all(results.repetition < 0.5, :) =  - vbias_all(results.repetition < 0.5, :);
    boundedline(cohlevels, nanmean(vbias_all), nanstd(vbias_all) ./ sqrt(length(results.repetition)));
    set(gca, 'xtick', cohlevels, 'xticklabelrotation', -30);
    xlabel('Coherence (%)');
    ylabel({'Drift bias' '(in direction of p(repeat)'});
    ylim([0 0.4]);
    box off; offsetAxes;
    title(datasetnames{d});
    
    
end


tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/vbias_per_coherence.pdf'));
