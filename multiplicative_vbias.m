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
