function alldat = individual_correlation_prevcorrect
% from the csv table, make an overview of repetition behaviour

% get a huge list with values for each participant
% can then work with this dataframe

close all; clc;
addpath(genpath('~/code/Tools'));

global mypath datasets datasetnames colors
cnt = 1;

if ~exist('Gsq', 'var'), Gsq = 0; end
if ~exist('sz', 'var'),  sz = 0; end

% ============================================ %
% ONE LARGE PLOT WITH PANEL FOR EACH DATASET
% ============================================ %

doText = false;

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
    % ============================================ %
    
    % use the stimcoding difference
    results.z_prevresp_correct = ...
        results.(['z_1_1__' whichmdls 'dczprevcorrect']) - results.(['z_1_2__' whichmdls 'dczprevcorrect']);
    results.v_prevresp_correct = ...
        results.(['dc_1_1__' whichmdls 'dczprevcorrect']) - results.(['dc_1_2__' whichmdls 'dczprevcorrect']);
    
    try
        results.z_prevresp_error = ...
            results.(['z_0_1__' whichmdls 'dczprevcorrect']) - results.(['z_0_2__' whichmdls 'dczprevcorrect']);
        results.v_prevresp_error = ...
            results.(['dc_0_1__' whichmdls 'dczprevcorrect']) - results.(['dc_0_2__' whichmdls 'dczprevcorrect']);
    catch
        results.z_prevresp_error = ...
            results.z_c10__stimcodingdczprevcorrect - results.(['z_0_2__' whichmdls 'dczprevcorrect']);
        results.v_prevresp_error = ...
            results.dc_c10__stimcodingdczprevcorrect - results.(['dc_0_2__' whichmdls 'dczprevcorrect']);
    end
        
    cols = cbrewer('qual', 'Paired', 10);
    
    % assign to structure
    allresults(1).z_prevresp     = results.z_prevresp_correct;
    allresults(1).v_prevresp     = results.v_prevresp_correct;
    allresults(1).criterionshift = results.repetition;
    alltitles{1}                 = cat(2, datasetnames{d}{1}, ' - ', 'Correct');
    allresults(1).marker 			= 'o';
    allresults(1).meancolor 		= [ 0 0 0];
    allresults(1).scattercolor	 	= [ 0.5 0.5 0.5];
    
    % also after error choices
    allresults(2).z_prevresp     = results.z_prevresp_error;
    allresults(2).v_prevresp     = results.v_prevresp_error;
    allresults(2).criterionshift = results.repetition;
    alltitles{2}                 = cat(2, datasetnames{d}{1}, ' - ', 'Error');
    allresults(2).marker        = 's';
    allresults(2).meancolor 		= cols(6, :);
    allresults(2).scattercolor	 	= cols(5, :);
    
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
        'rows', 'complete', 'type', 'spearman');
    if pval3 < 0.05,
        fprintf('warning %s: rho = %.3f, pval = %.3f \n', datasets{d}, rho3, pval3);
    end
    [rhodiff, ~, pval] = rddiffci(rho1,rho2,rho3,numel(~isnan( cat(1, allresults(:).criterionshift))), 0.05);
    
    % move together
    sp2.Position(1) = sp2.Position(1) - 0.08;
    ss = suplabel(cat(2, datasetnames{d}{1}, ' ', datasetnames{d}{2}), 't');
    set(ss, 'fontweight', 'normal');
    ss.FontWeight = 'normal';
    ss.Position(2) = ss.Position(2) - 0.03;
    
    % add colored axes after suplabel (which makes them black)
    xlabel(sp1, 'History shift in z');
    set(sp1, 'xcolor', colors(1, :));
    xlabel(sp2, 'History shift in v_{bias}');
    set(sp2, 'xcolor', colors(2, :));
    
    if doText,
        %% add line between the two correlation coefficients
        txt = {sprintf('\\Delta\\rho(%d) = %.3f, p = %.3f', length(find(~isnan(cat(1, allresults(:).criterionshift) )))-3, rhodiff, pval)};
        if pval < 0.001,
            txt = {sprintf('\\Delta\\rho(%d) = %.3f, p < 0.001', length(find(~isnan(cat(1, allresults(:).criterionshift) )))-3,  rhodiff)};
        end
        title(txt, 'fontweight', 'normal', 'fontsize', 8, 'horizontalalignment', 'left');
    end
    
    tightfig;
    if Gsq,
        print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/figure1c_Gsq_modelfree_prevcorrect_sz%d_d%d.pdf', d, sz));
    else
        print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/figure1c_HDDM_modelfree_prevcorrect_sz%d_d%d.pdf', d, sz));
        %print(gcf, '-depsc', sprintf('~/Data/serialHDDM/figure1c_HDDM_modelfree_stimcoding_d%d.eps', d));
    end
    
    for a = 1:length(allresults),
        
        % % SAVE CORRELATIONS FOR OVERVIEW PLOT
        % [r,p,rlo,rup] = corrcoef(allresults(a).z_prevresp, allresults(a).criterionshift);
        % alldat(cnt).corrz = r(1,2);
        % alldat(cnt).corrz_ci = [rlo(1,2) rup(1,2)];
        % alldat(cnt).pz = p(1,2);
        % alldat(cnt).bfz = corrbf(r(1,2), numel(allresults(a).z_prevresp));
        
        % [r,p,rlo,rup] = corrcoef(allresults(a).v_prevresp, allresults(a).criterionshift);
        % alldat(cnt).corrv = r(1,2);
        % alldat(cnt).corrv_ci = [rlo(1,2) rup(1,2)];
        % alldat(cnt).pv = p(1,2);
        % alldat(cnt).bfv = corrbf(r(1,2), numel(allresults(a).v_prevresp));
        

           % SAVE CORRELATIONS FOR OVERVIEW PLOT
        % COMPUTE THE SPEARMANS CORRELATION AND ITS CONFIDENCE INTERVAL!
        [alldat(cnt).corrz, alldat(cnt).corrz_ci, alldat(cnt).pz, alldat(cnt).bfz] = ...
            spearmans(allresults(a).z_prevresp, allresults(a).criterionshift);

        [alldat(cnt).corrv, alldat(cnt).corrv_ci, alldat(cnt).pv, alldat(cnt).bfv] = ...
            spearmans(allresults(a).v_prevresp, allresults(a).criterionshift);

        alldat(cnt).datasets = datasets{d};
        alldat(cnt).datasetnames = datasetnames{d};
        
        % also add the difference in r, Steigers test
        [r,p,rlo,rup] = corrcoef(allresults(a).v_prevresp, allresults(a).z_prevresp);
        
        [rhodiff, rhodiffci, pval] = rddiffci(alldat(cnt).corrz,alldat(cnt).corrv, ...
            r(1,2), ...
            numel(allresults(a).v_prevresp), 0.05);
        
        alldat(cnt).corrdiff = rhodiff;
        alldat(cnt).corrdiff_ci = rhodiffci;
        alldat(cnt).pdiff = pval;
        alldat(cnt).nsubj  = numel(allresults(a).criterionshift);
        
            
        alldat(cnt).marker          = allresults(a).marker;
        alldat(cnt).scattercolor    = allresults(a).scattercolor;
        alldat(cnt).meancolor       = allresults(a).meancolor;
        
        cnt = cnt + 1;
    end
end

end
