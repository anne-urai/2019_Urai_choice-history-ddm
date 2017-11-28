function e6_serialBias_SfN_modelFree_CRF_PPC
% ========================================== %
% conditional response functions from White & Poldrack
% run on simulated rather than real data
% ========================================== %

addpath(genpath('~/code/Tools'));
warning off; close all; clear;
global datasets datasetnames mypath colors

qntls{1} = [.2, .4, .6, .8, .95]; % White & Poldrack
% qntls = [.1, .3, .5, .7, .9, 1]; % Leite & Ratcliff
qntls{2} = [0.1, 0.3, 0.5, 0.7, 0.9];

allcols = colors;

for d = 6:length(datasets),
	for q = 1:length(qntls),
	
    switch datasets{d}
        case {'Bharath_fMRI', 'Anke_MEG', 'Anke_2afc_sequential', 'Anke_merged'}
            tps = [20 80];
        otherwise
            tps = 0;
    end
    
    % plot
    close all;
    subplot(441); hold on;
    
    for tp = 1:length(tps),
        
        % redo this for each simulation
        models = {'stimcoding_z_prevresp', 'stimcoding_dc_prevresp', 'stimcoding_dc_z_prevresp', 'stimcoding_nohist', 'stimcoding_nohist'};
        thesecolors = {colors(1, :), colors(2, :), {colors(1,:), colors(2, :)}, [0.5 0.5 0.5], [0 0 0]};
        
        for m = 1:length(models),
            
            if ~exist(sprintf('%s/%s/%s/ppc_data.csv', mypath, datasets{d}, models{m}), 'file'),
                continue;
            else
                fprintf('%s/%s/%s/ppc_data.csv \n', mypath, datasets{d}, models{m});
            end
            
            % load simulated data - make sure this has all the info we need
            alldata    = readtable(sprintf('%s/%s/%s/ppc_data.csv', mypath, datasets{d}, models{m}));
            alldata    = sortrows(alldata, {'subj_idx'});
            
            if ~any(ismember(alldata.Properties.VariableNames, 'transitionprob'))
                alldata.transitionprob = zeros(size(alldata.subj_idx));
            else
                assert(nanmean(unique(alldata.transitionprob)) == 50, 'rescale units');
                alldata = alldata(alldata.transitionprob == tps(tp), :);
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
					fprintf('splitting by coherence, nsubj %d newNsubj %d', origNsubj, newNsubj);
				end
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
            if tps(tp) == 0,
                alldata.biased(altIdx) = double(~(alldata.biased(altIdx))); % flip
		    end
            
            % divide RT into quantiles for each subject
            discretizeRTs = @(x) {discretize(x, quantile(x, [0, qntls{q}]))};
            rtbins = splitapply(discretizeRTs, alldata.rt, findgroups(alldata.subj_idx));
            alldata.rtbins = cat(1, rtbins{:});
            
            % get RT quantiles for choices that are in line with or against the bias
            [gr, sjidx, rtbins] = findgroups(alldata.subj_idx, alldata.rtbins);
            cpres               = array2table([sjidx, rtbins], 'variablenames', {'subj_idx', 'rtbin'});
            cpres.choice        = splitapply(@nanmean, alldata.biased, gr); % choice proportion
            
            % make into a subjects by rtbin matrix
            mat = unstack(cpres, 'choice', 'rtbin');
            mat = mat{:, 2:end}; % remove the last one, only has some weird tail
            
            % biased choice proportion
            if m < length(models),
                if isnumeric(thesecolors{m})
                    plot(qntls{q}, nanmean(mat, 1), 'color', thesecolors{m}, 'linewidth', 1);
                elseif iscell(thesecolors{m}) % superimposed lines for dashed
                    plot(qntls{q}, nanmean(mat, 1), 'color', thesecolors{m}{1}, 'linewidth', 1);
                    plot(qntls{q}, nanmean(mat, 1), ':', 'color', thesecolors{m}{2}, 'linewidth', 1);
                end
            else
                %% ALSO ADD THE REAL DATA
                h = ploterr(qntls{q}, nanmean(mat, 1), [], nanstd(mat, [], 1) ./ sqrt(size(mat, 1)), 'k-', 'abshhxy', 0);
                set(h(1), 'color', 'k', 'marker', '.', 'markerfacecolor', 'k', 'markeredgecolor', 'k', 'linewidth', 0.5, 'markersize', 10);
                set(h(2), 'linewidth', 0.5);
            end
        end
	end
		
        axis tight; box off;
        set(gca, 'xtick', qntls{q});
        axis square;  offsetAxes;
        xlabel('RT (quantiles)');
        
        if tps(tp) > 0,
            switch tps(tp)
                case 20
                    title([datasetnames{d}{1} ' Alternating']);
                case 50
                    title([datasetnames{d}{1} ' Neutral']);
                case 80
                    title([datasetnames{d}{1} ' Repetitive']);
            end
            ylabel('Fraction repetitions');
        else
            ylabel('Fraction biased choices');
        end
        
        title(datasetnames{d}{1});
        tightfig;
        print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/CRF_PPC_d%d_qntlsR%d.pdf', d, q));
        fprintf('~/Data/serialHDDM/CRF_PPC_d%d.pdf \n', d);
    
end
end
end
