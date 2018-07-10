function e6_serialBias_CBFs_perCoherence

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
models = {'data', 'stimcoding_nohist', 'stimcoding_z_prevresp',  ...
    'stimcoding_dc_prevresp', 'stimcoding_dc_z_prevresp'};
thesecolors = {[0 0 0], [0.5 0.5 0.5],  colors(1, :), ...
    colors(2, :), mean(colors([1 2], :))};

% only use the 2 datasets with coherence levels
for d = [4 5];
    
    % plot
    close all;
    figure; hold on;
    
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
                
                alldata.rt          = abs(alldata.rt_sampled);
                alldata.response    = alldata.response_sampled;
        end
        
        % make sure to use absolute RTs!
        alldata.rt = abs(alldata.rt);
        
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
        
        % divide RT into quantiles for each subject
        discretizeRTs = @(x) {discretize(x, quantile(x, [0, qntls]))};
        alldata(isnan(alldata.rt), :) = [];
        
        disp('splitting by coherence first');
        cohs = unique(alldata.coherence);
        for c = 1:length(cohs),
            
            subplot(3,3,c); hold on;
            thisdata = alldata(alldata.coherence == cohs(c), :);
            
            [gr, sj] = findgroups(thisdata.subj_idx);
            rtbins = splitapply(discretizeRTs, thisdata.rt, gr);
            thisdata.rtbins = cat(1, rtbins{:});
            
            % get RT quantiles for choices that are in line with or against the bias
            [gr, sjidx, rtbins, coh] = findgroups(thisdata.subj_idx, thisdata.rtbins, thisdata.coherence);
            cpres               = array2table([sjidx, rtbins], 'variablenames', {'subj_idx', 'rtbin'});
            cpres.choice        = splitapply(@nanmean, thisdata.biased, gr); % choice proportion
           
            % make into a subjects by rtbin matrix
            mat = unstack(cpres, 'choice', 'rtbin');
            mat = mat{:, 2:end}; % remove the last one, only has some weird tail
            
            % biased choice proportion
            switch models{m}
                case 'data'
                    % ALSO ADD THE REAL DATA WITH SEM
                    h = ploterr(qntls, nanmean(mat, 1), [], ...
                        nanstd(mat, [], 1) ./ sqrt(size(mat, 1)), 'k', 'abshhxy', 0);
                    set(h(1), 'color', 'k', 'marker', '.', ...
                        'markerfacecolor', 'k', 'markeredgecolor', 'k', 'linewidth', 0.5, 'markersize', 10, ...
                        'linestyle', '-');
                    set(h(2), 'linewidth', 0.5);
                otherwise
                    plot(qntls, nanmean(mat, 1), 'color', thesecolors{m}, 'linewidth', 1);
            end
            
            axis tight; box off;
            set(gca, 'xtick', qntls);
            axis square;  offsetAxes;
            xlabel('RT (quantiles)');
            set(gca, 'xcolor', 'k', 'ycolor', 'k');
            ylabel('P(bias)');
            set(gca, 'xcolor', 'k', 'ycolor', 'k');
            if rem(coh(c) * 100, 1) == 0,
                title(sprintf('Coherence %d%%', coh(c) * 100));
            else
                title(sprintf('Coherence %.3f%%', coh(c) * 100));
            end
        end
    end
    
    suplabel(cat(2, datasetnames{d}{1}, ' ', datasetnames{d}{2}), 't');
    tightfig;
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/CBF_d%d_perCoherence.pdf', d));
end