function conditional_bias_functions_collapsed(whichModels, qidx, xAxis, useBiasedSj, subject_cutoff)

% Code to fit the history-dependent drift diffusion models described in
% Urai AE, Gee JW de, Donner TH (2018) Choice history biases subsequent evidence accumulation. bioRxiv:251595
%
% MIT License
% Copyright (c) Anne Urai, 2018
% anne.urai@gmail.com

% ========================================== %
% conditional response functions from White & Poldrack
% run on simulated rather than real data
% ========================================== %

addpath(genpath('~/code/Tools'));
warning off; % close all;
global datasets datasetnames mypath colors

groups = {'alternators', 'repeaters', 'all'};
g = 3; % combine all subjects

switch qidx
    case 1
        qntls = [0.1, 0.3, 0.5, 0.7, 0.9, 1]; % Leite & Ratcliff
    case 2
        qntls = [0.5 1];
    case 3
        qntls = [.2, .4, .6, .8, .95]; % white & poldrack
end

switch whichModels
    case 1
        models = {'data', 'stimcoding_nohist', 'stimcoding_dc_z_prevresp'};
        thesecolors = {[0 0 0], [0.5 0.5 0.5], mean(colors([1 2], :))};
    case 2
        models = {'data', 'stimcoding_z_prevresp', 'stimcoding_dc_z_prevresp'};
        thesecolors = {[0 0 0],  colors(1, :), mean(colors([1 2], :))};
    case 3
        models = {'data', 'stimcoding_z_prevresp', 'stimcoding_dc_prevresp', 'stimcoding_dc_z_prevresp'};
        thesecolors = {[0 0 0], colors(1, :), colors(2, :), mean(colors([1 2], :))};
    case 4
        models = {'data', 'stimcoding_z_prevresp', 'stimcoding_dc_prevresp'};
        thesecolors = {[0 0 0], colors(1, :), colors(2, :)};
end

% ========================================== %
% START
% ========================================== %

% close all;
% subplot(3,3,1); 
hold on;

for m = 1:length(models),
    
    filename_save = sprintf('%s/summary/%s_ppc.mat', mypath, models{m});
    if exist(filename_save, 'file'),
        load(filename_save);
    else
        for d = 1:length(datasets)
            
            switch models{m}
                case 'data'
                    filename = dir(sprintf('%s/%s/*.csv', mypath, datasets{d}));
                    alldata_thisdata  = readtable(sprintf('%s/%s/%s', mypath, datasets{d}, filename.name));
                    
                otherwise
                    if ~exist(sprintf('%s/summary/%s/%s_ppc_data.csv', mypath, datasets{d}, models{m}), 'file'),
                        continue;
                    else
                        fprintf('%s/summary/%s/%s_ppc_data.csv \n', mypath, datasets{d}, models{m});
                    end
                    % load simulated data - make sure this has all the info we need
                    alldata_thisdata    = readtable(sprintf('%s/summary/%s/%s_ppc_data.csv', mypath, datasets{d}, models{m}));
                    alldata_thisdata    = sortrows(alldata_thisdata, {'subj_idx'});
                    
                    alldata_thisdata.rt          = abs(alldata_thisdata.rt_sampled);
                    alldata_thisdata.response    = alldata_thisdata.response_sampled;
            end
            
            % for those datasets with varying coherence, take only the difficult trials
            % FOR THE TWO DATASETS WITH VARYING COHERENCE, REMOVE THE HIGH LEVELS
            if any(ismember('coherence', alldata_thisdata.Properties.VariableNames)),
                if max(alldata_thisdata.coherence) == 81,
                    alldata_thisdata = alldata_thisdata(alldata_thisdata.coherence < 27, :);
                elseif max(alldata_thisdata.coherence) == 0.3,
                    alldata_thisdata = alldata_thisdata(alldata_thisdata.coherence < 0.1, :);
                else % do nothing
                end
            end
            
            % add some more stuff
            alldata_thisdata.subj_idx  = alldata_thisdata.subj_idx + 1000*d; % unique identifiers
            alldata_thisdata = alldata_thisdata(:, {'rt', 'response', 'prevresp', 'subj_idx'});
            alldata_thisdata.dataset(:, 1) = d;
            alldata_append{d} = alldata_thisdata;
        end
        alldata = cat(1, alldata_append{:});
        savefast(filename_save, 'alldata');
    end
    
    %  if m == 1, % preallocate
    %   allds.fast  = nan(numel(datasets), length(models));
    % 	allds.slow  = nan(numel(datasets), length(models));
    % end
    
    % make sure to use absolute RTs!
    alldata.rt = abs(alldata.rt);
    
    % recode into repeat and alternate for the model
    alldata.repeat = zeros(size(alldata.response));
    alldata.repeat(alldata.response == (alldata.prevresp > 0)) = 1;
    alldata.biased = alldata.repeat;
    
    % compute individual bias based on the data only
    switch models{m}
        case 'data'
            % for each observers, compute their bias
            [gr, sjs] = findgroups(alldata.subj_idx);
            sjrep = splitapply(@nanmean, alldata.repeat, gr);
            % who are the repeating observers?
            alternators = sjs(sjrep < 0.5);
    end
    
    %% subselect some of the observers
    switch groups{g}
        case 'alternators'
            alldata(~ismember(alldata.subj_idx, alternators), :) = [];
        case 'repeaters'
            alldata(ismember(alldata.subj_idx, alternators), :) = [];
        otherwise
            % recode into biased and unbiased choices
            altIdx = ismember(alldata.subj_idx, alternators);
            alldata.biased(altIdx) = double(~(alldata.biased(altIdx))); % flip
    end
    
    % if all subjects fall into one category..
    % specifically, JW's 1st dataset doesn't have any alternators
    if isempty(alldata),
        continue;
    end
    
    % select only a subset of the observers
    if useBiasedSj ~= 0,
        
        switch models{m}
            case 'data'
                % for this plot, use only the most extremely biased observers (see email Tobi 23 August)
                [gr, sjs] = findgroups(alldata.subj_idx);
                sjbias = splitapply(@nanmean, alldata.biased, gr);
                % take percentile
                cutoff = prctile(sjbias, 100-subject_cutoff);
                
                if useBiasedSj == 1,
                    usesj = sjs(sjbias > cutoff); % highest
                elseif useBiasedSj == -1,
                    usesj = sjs(sjbias < cutoff); % lowest
                end
        end
        alldata = alldata(ismember(alldata.subj_idx, usesj), :);
    end
    
    % divide RT into quantiles for each subject
    discretizeRTs = @(x) {discretize(x, quantile(x, [0, qntls]))};
    alldata(isnan(alldata.rt), :) = [];
    
    % discretize into bins of RT
    rtbins = splitapply(discretizeRTs, alldata.rt, findgroups(alldata.subj_idx));
    alldata.rtbins = cat(1, rtbins{:});
    
    % get RT quantiles for choices that are in line with or against the bias
    [gr, sjidx, rtbins] = findgroups(alldata.subj_idx, alldata.rtbins);
    cpres               = array2table([sjidx, rtbins], 'variablenames', {'subj_idx', 'rtbin'});
    cpres.choice        = splitapply(@nanmean, alldata.biased, gr); % choice proportion
    %cpres.meanrt        = splitapply(@nanmean, alldata.rt, gr); % average RT per bin
    
    % make into a subjects by rtbin matrix
    mat_tmp = unstack(cpres, 'choice', 'rtbin');
    mat     = mat_tmp{:, 2:end}; % remove the last one, only has some weird tail
    
    % also compute the mean RT for each subject and RT bin
    rtAvg               = array2table([sjidx, rtbins], 'variablenames', {'subj_idx', 'rtbin'});
    rtAvg.rt            = splitapply(@nanmean, alldata.rt, gr); % choice proportion
    xRTs                = unstack(rtAvg, 'rt', 'rtbin');
    xRTs                = xRTs{:, 2:end}; % remove the last one, only has some weird tail
    assert(isequal(size(mat), size(xRTs)), 'mismatch');

    % average within datasets first
    try
        mat  = splitapply(@nanmean, mat, findgroups(floor(mat_tmp.subj_idx / 1000)));
        % xRTs = splitapply(@nanmean, xRTs, findgroups(floor(mat_tmp.subj_idx / 1000)));
    catch
        assert(1==0)
    end

    % ================================ %
    % PLOT THE ACTUAL CURVES
    % ================================ %

    switch models{m}
        case 'data'
            switch xAxis
                case 'quantiles'
                    x_axis = qntls;
                    x_axis_std = [];
                case 'rt'
                    x_axis      = mean(xRTs);
                    ebar_ci     = bootci(5000, @mean, xRTs);
                    x_axis_std  = {ebar_ci(1, :), ebar_ci(2, :)}; % CI based on bootstrap
                    x_axis_std  = nanstd(xRTs) ./ sqrt(size(xRTs, 1)) * 1.96;
            end
            
            % bootstrapped CI
            ebar_ci = bootci(5000, @mean, mat);
            ebar_ci_sem = {ebar_ci(1, :), ebar_ci(2, :)}; % CI based on bootstrap
            ebar_ci_sem = nanstd(mat) ./ sqrt(size(mat, 1)) * 1.96; % CI based on s.e.m.
            
            % ALSO ADD THE REAL DATA WITH SEM/95%CI
            h = ploterr(x_axis, nanmean(mat, 1), x_axis_std, ebar_ci_sem, 'k', 'abshhxy', 0);
                        
            set(h(1), 'color', 'k', 'marker', 'o', ...
                'markerfacecolor', 'w', 'markeredgecolor', 'k', 'linewidth', 0.5, 'markersize', 3, ...
                'linestyle', '-');
            set([h(2) h(3)], 'linewidth', 0.5);
        otherwise
            plot(x_axis, nanmean(mat, 1), 'color', thesecolors{m}, 'linewidth', 1);
    end
    
    % SAVE lowest and highest quantiles
    allds.fast(:, m) = mat(:, 1);
    allds.slow(:, m) = mat(:, end);
end

axis tight; box off;
set(gca, 'xtick', roundn(x_axis, -2), 'xticklabelrotation', -30);
% ylim([0.5 0.56]);

axis square;  offsetAxes;
switch xAxis
    case 'quantiles'
        xlabel('Response time (quantile)')
    case 'rt'
        xlabel('Response time (s)');
end

% switch useBiasedSj
%     case 0
%         title('All subjects');
%     case -1
%         title(sprintf('Least biased %d percentile of subjects', 100-subject_cutoff));
%     case 1
%         title(sprintf('Most biased %d percentile of subjects', subject_cutoff));
% end

set(gca, 'xcolor', 'k', 'ycolor', 'k');
ylabel(sprintf('P(bias), %s', groups{g}));
ylabel('Choice bias (fraction)');
% title('Collapsed across datasets');
set(gca, 'xcolor', 'k', 'ycolor', 'k');
% 
% tightfig;
% 
% print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/CBFs_q%d_%s_sjCutoff%d_%dpercentile_models%d.pdf', ...
%     qidx, xAxis, useBiasedSj, subject_cutoff, whichModels));
% 
% % print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/CBFs_collapsed_q%d_%s_sjCutoff%d_%s_models%d.pdf', ...
% %  qidx, xAxis, useBiasedSj, groups{g}, whichModels));

%% ========================================== %
% PLOT ACROSS DATASETS - only for median split
% ========================================== %

savefast(sprintf('~/Data/serialHDDM/allds_cbfs_mediansplit.mat'), 'allds', 'models');
disp('saved output');


end