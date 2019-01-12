
function alldat = individual_correlation_tcoh()

% Code to fit the history-dependent drift diffusion models described in
% Urai AE, Gee JW de, Donner TH (2018) Choice history biases subsequent evidence accumulation. bioRxiv:251595
%
% MIT License
% Copyright (c) Anne Urai, 2018
% anne.urai@gmail.com

close all; clc;
addpath(genpath('~/code/Tools'));

global mypath datasets datasetnames colors
cnt = 1;

if ~exist('Gsq', 'var'), Gsq = 0; end
if ~exist('sz', 'var'),  sz = 0; end

% ============================================ %
% ONE LARGE PLOT WITH PANEL FOR EACH DATASET
% ============================================ %

doText = true;

switch sz
    case 1
        whichmdls = ['stimcodingst'];
    case 0
        whichmdls = ['stimcoding'];
end

% color in different grouos
tmpcolors = cbrewer('qual', 'Paired', 10);
transitioncolors = [[0.5 0.5 0.5]; tmpcolors([7 9], :)];
meancolors = [0 0 0; tmpcolors([8 10], :)];
markers = {'o', 'v', '^'}; %also indicate with different markers

close all;
for d = 2:3,
    disp(datasets{d});
    
    % colors = [8 141 165; 141 165 8;  150 150 150] ./ 256;
    
    if Gsq,
        results = readtable(sprintf('%s/summary/%s/allindividualresults_Gsq.csv', mypath, datasets{d}));
    else
        results = readtable(sprintf('%s/summary/%s/allindividualresults.csv', mypath, datasets{d}));
    end
    results = results(results.session == 0, :);
    
    allresults = struct(); alltitles = {};
    results.z_prevresp = ...
        results.z_1__stimcodingdczprevrespstcoh - results.z_2__stimcodingdczprevrespstcoh;
    results.v_prevresp = ...
        results.dc_1__stimcodingdczprevrespstcoh - results.dc_2__stimcodingdczprevrespstcoh;
    results.criterionshift = results.repetition;
    
    % assign to structure
    allresults(1).z_prevresp     = results.z_prevresp;
    allresults(1).v_prevresp     = results.v_prevresp;
    allresults(1).criterionshift = results.criterionshift;
    
    allresults(1).marker 			= markers{1};
    allresults(1).meancolor 		= meancolors(1, :);
    allresults(1).scattercolor	 	= transitioncolors(1, :);
    alltitles{1} 					= {datasetnames{d}{1} datasetnames{d}{2}}; % use only the dataset title
    

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
    try
        ss = suplabel(cat(2, datasetnames{d}{1}, ' ', datasetnames{d}{2}), 't');
    catch
        ss = suplabel(datasetnames{d}{1}, 't');
    end
    
    set(ss, 'fontweight', 'normal');
    ss.FontWeight = 'normal';
    ss.Position(2) = ss.Position(2) - 0.03;
    
    % add colored axes after suplabel (which makes them black)
    xlabel(sp1, 'History shift in z');
    set(sp1, 'xcolor', colors(1, :), 'ycolor', 'k');
    xlabel(sp2, 'History shift in v_{bias}');
    set(sp2, 'xcolor', colors(2, :), 'ycolor', 'k');
    
    if doText,
        %% add line between the two correlation coefficients
        txt = {sprintf('\\Delta\\rho(%d) = %.3f, p = %.3f', length(find(~isnan(cat(1, allresults(:).criterionshift) )))-3, rhodiff, pval)};
        if pval < 0.001,
            txt = {sprintf('\\Delta\\rho(%d) = %.3f, p < 0.001', length(find(~isnan(cat(1, allresults(:).criterionshift) )))-3,  rhodiff)};
        end
        tt = title(txt, 'fontweight', 'normal', 'fontsize', 6, 'horizontalalignment', 'left');
        tt.Position(2) = tt.Position(2) - 0.008;
    end
    
    tightfig;
    if Gsq,
        print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/figure1c_Gsq_modelfree_stimcoding_tcoh_d%d.pdf', d));
    else
        print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/figure1c_HDDM_modelfree_stimcoding_tcoh_d%d.pdf', d ));
        %print(gcf, '-depsc', sprintf('~/Data/serialHDDM/figure1c_HDDM_modelfree_stimcoding_d%d.eps', d));
    end
    
    for a = 1:length(allresults),
        
        % SAVE CORRELATIONS FOR OVERVIEW PLOT
        % COMPUTE THE SPEARMANS CORRELATION AND ITS CONFIDENCE INTERVAL!
        [alldat(cnt).corrz, alldat(cnt).corrz_ci, alldat(cnt).pz, alldat(cnt).bfz] = ...
            spearmans(allresults(a).z_prevresp, allresults(a).criterionshift);

        [alldat(cnt).corrv, alldat(cnt).corrv_ci, alldat(cnt).pv, alldat(cnt).bfv] = ...
            spearmans(allresults(a).v_prevresp, allresults(a).criterionshift);

        % alldat(cnt).corrz = r(1,2);
        % alldat(cnt).corrz_ci = [rlo(1,2) rup(1,2)];
        % alldat(cnt).pz = p(1,2);
        % alldat(cnt).bfz = corrbf(r(1,2), numel(allresults(a).z_prevresp));
        
        % [r,p,rlo,rup] = corrcoef(allresults(a).v_prevresp, allresults(a).criterionshift);
        % alldat(cnt).corrv = r(1,2);
        % alldat(cnt).corrv_ci = [rlo(1,2) rup(1,2)];
        % alldat(cnt).pv = p(1,2);
        % alldat(cnt).bfv = corrbf(r(1,2), numel(allresults(a).v_prevresp));
        
        alldat(cnt).datasets        = datasets{d};
        alldat(cnt).datasetnames    = alltitles{a};
        
        % also add the difference in correlation, steigers test
        [r,p,rlo,rup]               = spearmans(allresults(a).v_prevresp, allresults(a).z_prevresp);
        
        [rhodiff, rhodiffci, pval] = rddiffci(alldat(cnt).corrz, alldat(cnt).corrv, ...
            r, numel(allresults(a).v_prevresp), 0.05);
        
        alldat(cnt).corrdiff        = rhodiff;
        alldat(cnt).corrdiff_ci     = rhodiffci;
        alldat(cnt).pdiff           = pval;
        
        % plotting layout for forestPlot
        alldat(cnt).marker          = allresults(a).marker;
        alldat(cnt).scattercolor    = allresults(a).scattercolor;
        alldat(cnt).meancolor       = allresults(a).meancolor;
        
        cnt = cnt + 1;
    end
end

end
