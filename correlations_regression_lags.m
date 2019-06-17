function allresults = correlations_regression_lags

% Code to fit the history-dependent drift diffusion models as described in
% Urai AE, de Gee JW, Tsetsos K, Donner TH (2019) Choice history biases subsequent evidence accumulation. eLife, in press.
%
% MIT License
% Copyright (c) Anne Urai, 2019
% anne.urai@gmail.com

global mypath datasets datasetnames
addpath(genpath('~/code/Tools'));
warning off; close all;
cols = cbrewer('qual', 'Paired', 10);

numlags = 6;
vars = {'z_correct', 'z_error', 'v_correct', 'v_error', 'repeat_correct', 'repeat_error'};
cnt = 1;

for d = 1:length(datasets),

    dat = readtable(sprintf('%s/summary/%s/allindividualresults.csv', mypath, datasets{d}));
    dat = dat(dat.session == 0, :);

    for m = 1:length(vars),
        alldata.(vars{m}) = nan(numlags, size(dat, 1));
    end

    % ALL MODELS THAT WERE RAN
    mdls = {'regress_nohist', ...
        'regress_z_lag1', ...
        'regress_dc_lag1', ...
        'regress_dcz_lag1', ...
        'regress_z_lag2', ...
        'regress_dc_lag2', ...
        'regress_dcz_lag2', ...
        'regress_z_lag3', ...
        'regress_dc_lag3', ...
        'regress_dcz_lag3', ...
        'regress_z_lag4', ...
        'regress_dc_lag4', ...
        'regress_dcz_lag4', ...
        'regress_z_lag5', ...
        'regress_dc_lag5', ...
        'regress_dcz_lag5', ...
        'regress_z_lag6', ...
        'regress_dc_lag6', ...
        'regress_dcz_lag6'};
    
    % ============================= %
    % 1. DETERMINE THE BEST MODEL
    % ============================= %
    
    mdldic = nan(1, length(mdls));
    for m = 1:length(mdls),
    try
        modelcomp = readtable(sprintf('%s/%s/%s/model_comparison.csv', ...
            mypath, datasets{d}, mdls{m}), 'readrownames', true);
        mdldic(m) = modelcomp.aic;
    catch
        fprintf('%s/%s/%s/model_comparison.csv  NOT FOUND\n', ...
            mypath, datasets{d}, mdls{m})
        end
    end

    % everything relative to the full model
    mdldic = bsxfun(@minus, mdldic, mdldic(1));
    mdldic = mdldic(2:end);
    mdls = mdls(2:end);
    [~, bestMdl] = min(mdldic);

    % now take the hybrid model for this best-fitting lag
    bestmodelname = sprintf('regressdczlag%s', mdls{bestMdl}(end));
    disp(bestmodelname);

    % ========================================================== %
    % 2. FOR THE BEST-FITTING MODEL, GET HISTORY WEIGHTS
    % ========================================================== %

    % ignore lag 1 - just take the average of lag 2:bestmodel
    for l = 2:str2double(bestmodelname(end)),
        lname = num2str(l);
                
        % get regression weights
        for v = 1:length(vars),
            switch vars{v}
            case 'z_correct'
            alldata.(vars{v})(l,:) = ...
                (dat.(['z_prev' lname 'resp__' bestmodelname]) + ...
                dat.(['z_prev' lname  'stim__' bestmodelname]));
            case 'z_error'
            alldata.z_error(l,:) = ...
                (dat.(['z_prev' lname  'resp__' bestmodelname]) - ...
                dat.(['z_prev' lname  'stim__' bestmodelname]));
            case 'v_correct'
            alldata.v_correct(l,:) = ...
                (dat.(['v_prev' lname  'resp__' bestmodelname]) + ...
                dat.(['v_prev' lname  'stim__' bestmodelname]));
            case 'v_error'
            alldata.v_error(l,:) = ...
                (dat.(['v_prev' lname  'resp__' bestmodelname]) - ...
                dat.(['v_prev' lname  'stim__' bestmodelname]));

            case 'repeat_error'
                alldata.(vars{v})(l,:) =  dat.(['repetition_error' num2str(l)])...
                - arrayfun(@trivial_probabilities, dat.repetition_error1, repmat(l, size(dat, 1), 1));

            case 'repeat_correct'
                alldata.(vars{v})(l,:) =  dat.(['repetition_correct' num2str(l)]) ...
                - arrayfun(@trivial_probabilities, dat.repetition_correct1, repmat(l, size(dat, 1), 1));
            end

        end 

    end

    % assign to structure - correct choices
    allresults(1).z_prevresp      = nanmean(alldata.z_correct);
    allresults(1).v_prevresp      = nanmean(alldata.v_correct);
    allresults(1).criterionshift  = nanmean(alldata.repeat_correct);
    alltitles{1}                  = datasetnames{d};
    allresults(1).marker            = 'o';
    allresults(1).meancolor         = [0 0 0];
    allresults(1).scattercolor      = [0.5 0.5 0.5];

    % also after error choices
    allresults(2).z_prevresp     = nanmean(alldata.z_error);
    allresults(2).v_prevresp     = nanmean(alldata.v_error);
    allresults(2).criterionshift = nanmean(alldata.repeat_error);
    alltitles{2}                 = datasetnames{d};
    allresults(2).marker        = 's';
    allresults(2).meancolor         = cols(6, :);
    allresults(2).scattercolor      = cols(5, :);

    % ========================================================== %
    % COMPUTE CORRELATIONS
    % ========================================================== %

    for a = 1:length(allresults),
     
     % SAVE CORRELATIONS FOR OVERVIEW PLOT
     % COMPUTE THE SPEARMANS CORRELATION AND ITS CONFIDENCE INTERVAL!
     [alldat(cnt).corrz, alldat(cnt).corrz_ci, alldat(cnt).pz, alldat(cnt).bfz] = ...
         spearmans(allresults(a).z_prevresp(:), allresults(a).criterionshift(:));

     [alldat(cnt).corrv, alldat(cnt).corrv_ci, alldat(cnt).pv, alldat(cnt).bfv] = ...
         spearmans(allresults(a).v_prevresp(:), allresults(a).criterionshift(:));

     alldat(cnt).datasets        = datasets{d};
     alldat(cnt).datasetnames    = alltitles{a};
     
     % also add the difference in correlation, steigers test
     [r,p,rlo,rup]               = spearmans(allresults(a).v_prevresp(:), allresults(a).z_prevresp(:));
     
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

% ========================================================== %
% COMPUTE CORRELATIONS
% ========================================================== %

    forestPlot(alldat(1:2:end));
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/forestplot_regressionHDDM_prevcorrect.pdf'));
    forestPlot(alldat(2:2:end));
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/forestplot_regressionHDDM_preverror.pdf'));

end

function [ vec_repeat ] = trivial_probabilities(p,lag)
vec_repeat(1)=p;
for i=2:lag;
    vec_repeat(i)=p*vec_repeat(i-1)+(1-p)*(1-vec_repeat(i-1));
end

vec_repeat = vec_repeat(end);
end

