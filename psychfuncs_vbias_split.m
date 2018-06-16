function psychfuncs_vbias_split

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
warning off; close all; clear;
global datasets datasetnames mypath colors
qntls = [0.1, 0.3, 0.5, 0.7, 0.9]; % Leite & Ratcliff

% redo this for each simulation
models = {'stimcoding_nohist', 'stimcoding_z_prevresp',  ...
    'stimcoding_dc_prevresp', 'stimcoding_dc_z_prevresp' ...
    'data'};
thesecolors = {[0.5 0.5 0.5],  colors(1, :), ...
    colors(2, :), mean(colors([1 2], :)), [0 0 0]};

fixedEffects = 0;
allcols = colors;

allds.fast  = nan(length(datasets), length(models));
allds.slow = nan(length(datasets), length(models));

for d = 1:length(datasets); %:-1:1,
    
    % plot
    close all;
    subplot(441); hold on;
    
    for m = 1:length(models),
        
        switch models{m}
            case 'data'
                filename = dir(sprintf('%s/%s/*.csv', mypath, datasets{d}));
                alldata  = readtable(sprintf('%s/%s/%s', mypath, datasets{d}, filename.name));
            otherwise
                if ~exist(sprintf('%s/summary/%s/%s_ppc_data.csv', mypath, datasets{d}, models{m}), 'file'),
                    continue;
                else
                    fprintf('%s/summary/%s/%s_ppc_data.csv \n', mypath, datasets{d}, models{m});
                end
                % load simulated data - make sure this has all the info we need
                alldata    = readtable(sprintf('%s/summary/%s/%s_ppc_data.csv', mypath, datasets{d}, models{m}));
                alldata    = sortrows(alldata, {'subj_idx'});
        end
        
        if m < length(models),
            % use the simulations rather than the subjects' actual responses
            alldata.rt          = abs(alldata.rt_sampled);
            alldata.response    = alldata.response_sampled;
        end
        
        % when there were multiple levels of evidence, do these plots
        % separately for each level
        if any(ismember(alldata.Properties.VariableNames, 'coherence'))
            origNsubj = numel(unique(alldata.subj_idx));
            alldata.subj_idx = findgroups(alldata.subj_idx, alldata.coherence);
            newNsubj = numel(unique(alldata.subj_idx));
            if origNsubj ~= newNsubj,
                fprintf('splitting by coherence, nsubj %d newNsubj %d \n', origNsubj, newNsubj);
            end
        end
        
        % make sure to use absolute RTs!
        alldata.rt = abs(alldata.rt);
        alldata.correct = ((alldata.stimulus > 0) == alldata.response);
        assert(nanmean(alldata.correct) > 0.5);
        
        % recode into repeat and alternate for the model
        alldata.repeat = zeros(size(alldata.response));
        alldata.repeat(alldata.response == (alldata.prevresp > 0)) = 1;
        
        % for each observers, compute their bias
        [gr, sjs] = findgroups(alldata.subj_idx);
        sjrep = splitapply(@nanmean, alldata.repeat, gr);
        sjrep = sjs(sjrep < 0.5);
        
        % recode into biased and unbiased choices
        alldata.biased = alldata.repeat;
        altIdx = ismember(alldata.subj_idx, sjrep);
        alldata.biased(altIdx) = double(~(alldata.biased(altIdx))); % flip
        
        % fixed effects
        if fixedEffects,
            alldata.subj_idx = ones(size(alldata.subj_idx));
        end
        
        % divide RT into quantiles for each subject
        discretizeRTs = @(x) {discretize(x, quantile(x, [0, qntls]))};
        rtbins = splitapply(discretizeRTs, alldata.rt, findgroups(alldata.subj_idx));
        alldata.rtbins = cat(1, rtbins{:});
        
        % get RT quantiles for choices that are in line with or against the bias
        [gr, sjidx, rtbins] = findgroups(alldata.subj_idx, alldata.rtbins);
        cpres               = array2table([sjidx, rtbins], 'variablenames', {'subj_idx', 'rtbin'});
        cpres.choice        = splitapply(@nanmean, alldata.correct, gr); % choice proportion
        
        % make into a subjects by rtbin matrix
        mat = unstack(cpres, 'choice', 'rtbin');
        mat = mat{:, 2:end}; % remove the last one, only has some weird tail
        
        % biased choice proportion
        if m < length(models),
            plot(qntls, nanmean(mat, 1), 'color', thesecolors{m}, 'linewidth', 1);
        else
            %% ALSO ADD THE REAL DATA
            h = ploterr(qntls, nanmean(mat, 1), [], ...
                1.96 *  nanstd(mat, [], 1) ./ sqrt(size(mat, 1)), 'k', 'abshhxy', 0);
            set(h(1), 'color', 'k', 'marker', '.', ...
                'markerfacecolor', 'k', 'markeredgecolor', 'k', 'linewidth', 0.5, 'markersize', 10, ...
                'linestyle', '-');
            set(h(2), 'linewidth', 0.5);
        end
        
        % SAVE
        avg = nanmean(mat, 1);
        allds.fast(d, m) = nanmean(avg(1:2));
        allds.slow(d, m) = nanmean(avg(end-3:end));
        allds.all(d, m, :) = avg;
    end
    %  end
    
    axis tight; box off;
    set(gca, 'xtick', qntls);
    axis square;  offsetAxes;
    xlabel('RT (quantiles)');
    set(gca, 'xcolor', 'k', 'ycolor', 'k');
    ylabel('P(correct)');
    
    title(datasetnames{d});
    tightfig;
    set(gca, 'xcolor', 'k', 'ycolor', 'k');
    
    if fixedEffects,
        print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/CAF_PPC_d%d_fixed.pdf', d));
    else
        print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/CAF_PPC_d%d_q2.pdf', d));
    end
    fprintf('~/Data/serialHDDM/CRF_PPC_d%d.pdf \n', d);
end

end
