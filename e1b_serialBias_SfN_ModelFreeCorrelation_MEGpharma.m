function alldat = e1b_serialBias_SfN_ModelFreeCorrelation_MEGpharma
% from the csv table, make an overview of repetition behaviour

close all; clc;
addpath(genpath('~/code/Tools'));

global mypath datasets datasetnames colors

datasets = {'MEG'}; % CHANGE THIS TO MEG_750MS ONCE FINISHED!
% datasets = {'MEG_MEGsessions'}; % CHANGE THIS TO MEG_750MS ONCE FINISHED!
datasetnames = {{'Visual motion' '2IFC (FD) #2'}};
cnt = 1;

Gsq = 0; sz = 0;
doText = true;

% ============================================ %
% ONE LARGE PLOT WITH PANEL FOR EACH DATASET
% ============================================ %

switch sz
    case 1
        whichmdls = ['stimcodingsz'];
    case 0
        whichmdls = ['stimcoding'];
end

close all;
for d = length(datasets):-1:1
    disp(datasets{d});
        
    if Gsq,
        results = readtable(sprintf('%s/summary/%s/allindividualresults_Gsq.csv', mypath, datasets{d}));
    else
        results = readtable(sprintf('%s/summary/%s/allindividualresults.csv', mypath, datasets{d}));
    end
    results = results(results.session == 0, :);
    allresults = struct(); alltitles = {};
    
    % ============================================ %
    % RECODE INTO HISTORY SHIFT POINT ESTIMATES
    % SEPARATELY FOR EACH PHARMA GROUP!
    % ============================================ %
    
    % color in different groups
    colors1 = cbrewer('qual', 'Set2', 8);
    colors2 = cbrewer('qual', 'Dark2', 8);
    paired = cbrewer('qual', 'Paired', 12);

    % use the stimcoding difference
    results.z_prevresp_placebo = ...
        results.(['z_placebo_1__' whichmdls 'dczprevresppharma']) - results.(['z_placebo_2__' whichmdls 'dczprevresppharma']);
    results.v_prevresp_placebo = ...
        results.(['dc_placebo_1__' whichmdls 'dczprevresppharma']) - results.(['dc_placebo_2__' whichmdls 'dczprevresppharma']);
    
    allresults(1).z_prevresp     = results.z_prevresp_placebo;
    allresults(1).v_prevresp     = results.v_prevresp_placebo;
    allresults(1).criterionshift = results.repetition;
    alltitles{1}                 = {cat(2, datasetnames{1}{1}, ' ', datasetnames{1}{2}) ' - Placebo'};
    allresults(1).scattercolor    = colors1(8, :);
    allresults(1).meancolor       = colors2(8, :);
    allresults(1).marker          = 'o';

    % ATOMOXETINE
      results.z_prevresp_atomoxetine = ...
        results.(['z_atomoxetine_1__' whichmdls 'dczprevresppharma']) - results.(['z_atomoxetine_2__' whichmdls 'dczprevresppharma']);
    results.v_prevresp_atomoxetine = ...
        results.(['dc_atomoxetine_1__' whichmdls 'dczprevresppharma']) - results.(['dc_atomoxetine_2__' whichmdls 'dczprevresppharma']);
    
    allresults(2).z_prevresp     = results.z_prevresp_atomoxetine;
    allresults(2).v_prevresp     = results.v_prevresp_atomoxetine;
    allresults(2).criterionshift = results.repetition;
    alltitles{2}                 = {cat(2, datasetnames{1}{1}, ' ', datasetnames{1}{2}) ' - Atomoxetine'};
    allresults(2).scattercolor    = paired(5, :);
    allresults(2).meancolor       = paired(6, :);
    allresults(2).marker          = 's';
    
    % DONEPEZIL
    results.z_prevresp_donepezil = ...
        results.(['z_donepezil_1__' whichmdls 'dczprevresppharma']) - results.(['z_donepezil_2__' whichmdls 'dczprevresppharma']);
    results.v_prevresp_donepezil = ...
        results.(['dc_donepezil_1__' whichmdls 'dczprevresppharma']) - results.(['dc_donepezil_2__' whichmdls 'dczprevresppharma']);
    
    allresults(3).z_prevresp     = results.z_prevresp_donepezil;
    allresults(3).v_prevresp     = results.v_prevresp_donepezil;
    allresults(3).criterionshift = results.repetition;
    alltitles{3}                 = {cat(2, datasetnames{1}{1}, ' ', datasetnames{1}{2}) ' - Donepezil'};
    allresults(3).scattercolor    = paired(1, :);
    allresults(3).meancolor       = paired(2, :);
    allresults(3).marker          = 'd';
    
    disp(datasets{d}); disp(numel(unique(results.subjnr)));
    close all;
    
    % PLOT
    sp1 = subplot(4,4,1); hold on;
    [rho1, tt1] = plotScatter(allresults, 'z_prevresp', 0.585, doText);
    ylabel('P(repeat)');
    
    sp2 = subplot(4,4,2); hold on;
    [rho2, tt2, handles] = plotScatter(allresults, 'v_prevresp', 0.05, doText);
    set(gca, 'yticklabel', []);
    
    set(sp2, 'ylim', get(sp1, 'ylim'), 'ytick', get(sp1, 'ytick'));
    
    % compute the difference in correlation
    [rho3, pval3] = corr(cat(1, allresults(:).v_prevresp), cat(1, allresults(:).z_prevresp), ...
        'rows', 'complete', 'type', 'pearson');
    if pval3 < 0.05,
        fprintf('warning %s: rho = %.3f, pval = %.3f \n', datasets{d}, rho3, pval3);
    end
    [rhodiff, ~, pval] = rddiffci(rho1,rho2,rho3,numel(~isnan( cat(1, allresults(:).criterionshift))), 0.05);
    
    % move together
    sp2.Position(1) = sp2.Position(1) - 0.08;
    ss = suplabel('Visual motion 2IFC (FD) #2', 't');
    set(ss, 'fontweight', 'normal');
    ss.FontWeight = 'normal';
    ss.Position(2) = ss.Position(2) - 0.03;
    
    % add colored axes after suplabel (which makes them black)
    xlabel(sp1, 'History bias in z');
    set(sp1, 'xcolor', colors(1, :));
    xlabel(sp2, 'History bias in v');
    set(sp2, 'xcolor', colors(2, :));
    
    if doText,
        %% add line between the two correlation coefficients
        txt = {sprintf('\\Deltar(%d) = %.3f, p = %.3f', length(find(~isnan(cat(1, allresults(:).criterionshift) )))-3, rhodiff, pval)};
        if pval < 0.001,
            txt = {sprintf('\\Deltar(%d) = %.3f, p < 0.001', length(find(~isnan(cat(1, allresults(:).criterionshift) )))-3,  rhodiff)};
        end
        title(txt, 'fontweight', 'normal', 'fontsize', 6, 'horizontalalignment', 'left');
    end

    tightfig;
    if Gsq,
        print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/figure1c_Gsq_modelfree_pharma_sz%d_d%d.pdf', d, sz));
    else
        print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/figure1c_HDDM_modelfree_pharma_sz%d_d%d.pdf', d, sz));
        %print(gcf, '-depsc', sprintf('~/Data/serialHDDM/figure1c_HDDM_modelfree_stimcoding_d%d.eps', d));
    end
    
    for a = 1:length(allresults),
        
        % SAVE CORRELATIONS FOR OVERVIEW PLOT
        [r,p,rlo,rup] = corrcoef(allresults(a).z_prevresp, allresults(a).criterionshift, 'rows', 'complete');
        alldat(cnt).corrz = r(1,2);
        alldat(cnt).corrz_ci = [rlo(1,2) rup(1,2)];
        alldat(cnt).pz = p(1,2);
        alldat(cnt).bfz = corrbf(r(1,2), numel(allresults(a).z_prevresp));
        
        [r,p,rlo,rup] = corrcoef(allresults(a).v_prevresp, allresults(a).criterionshift, 'rows', 'complete');
        alldat(cnt).corrv = r(1,2);
        alldat(cnt).corrv_ci = [rlo(1,2) rup(1,2)];
        alldat(cnt).pv = p(1,2);
        alldat(cnt).bfv = corrbf(r(1,2), numel(allresults(a).v_prevresp));
        
        alldat(cnt).datasets = datasets{d};
        alldat(cnt).datasetnames = alltitles{a};
        
        % also add the difference in r, Steigers test
        [r,p,rlo,rup] = corrcoef(allresults(a).v_prevresp, allresults(a).z_prevresp, 'rows', 'complete');
        
        [rhodiff, rhodiffci, pval] = rddiffci(alldat(cnt).corrz,alldat(cnt).corrv, ...
            r(1,2), numel(allresults(a).v_prevresp), 0.05);
        
        alldat(cnt).corrdiff = rhodiff;
        alldat(cnt).corrdiff_ci = rhodiffci;
        alldat(cnt).pdiff = pval;
        
        % some layout things
        alldat(cnt).meancolor = allresults(a).meancolor;
        alldat(cnt).scattercolor = allresults(a).scattercolor;
        alldat(cnt).marker = allresults(a).marker;
        
        cnt = cnt + 1;
    end
end

end
