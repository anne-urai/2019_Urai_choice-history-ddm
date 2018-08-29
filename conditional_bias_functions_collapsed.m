function conditional_bias_functions_collapsed(whichModels, qidx, useBiasedSj, g)

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
warning off; close all; 
global datasets datasetnames mypath colors
cutoff_quantile = 3;
groups = {'alternators', 'repeaters', 'all'};

switch qidx
case 1
    qntls = [0.1, 0.3, 0.5, 0.7, 0.9]; % Leite & Ratcliff
case 2
	qntls = [0.5 1];
end

switch whichModels
case 1
	models = {'data', 'stimcoding_nohist', 'stimcoding_dc_z_prevresp'};
	thesecolors = {[0 0 0], [0.5 0.5 0.5], mean(colors([1 2], :))};
case 2
	models = {'data', 'stimcoding_z_prevresp', 'stimcoding_dc_z_prevresp'};
	thesecolors = {[0 0 0],  colors(1, :), mean(colors([1 2], :))};
case 3
	models = {'data', 'stimcoding_nohist', 'stimcoding_z_prevresp', 'stimcoding_dc_prevresp', 'stimcoding_dc_z_prevresp'};
	thesecolors = {[0 0 0], [0.5 0.5 0.5],  colors(1, :), colors(2, :), mean(colors([1 2], :))};
end

% ========================================== %
% START
% ========================================== %

close all;
subplot(3,3,1); hold on;

        for m = 1:length(models),

        	for d = 1:length(datasets);

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

                % add some more stuff
                alldata_thisdata.subj_idx  = alldata_thisdata.subj_idx + 1000*d; % unique identifiers
                alldata_thisdata = alldata_thisdata(:, {'rt', 'response', 'prevresp', 'subj_idx'});
                alldata_thisdata.dataset(:, 1) = d;
                alldata_append{d} = alldata_thisdata;
            end
            alldata = cat(1, alldata_append{:});

            if m == 1, % preallocate
            	allds.fast  = nan(numel(datasets), length(models));
				allds.slow  = nan(numel(datasets), length(models));
			end

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
            if useBiasedSj,

            	switch models{m}
            	case 'data'
                % for this plot, use only the most extremely biased observers (see email Tobi 23 August)
                [gr, sjs] = findgroups(alldata.subj_idx);
                sjbias = splitapply(@nanmean, alldata.biased, gr);

                if cutoff_quantile == 2,
                	cutoff = median(sjbias);
                elseif cutoff_quantile == 0,
                        cutoff = 0; % keep everyone!
                    else
                    	cutoff = quantile(sjbias, cutoff_quantile);
                    end
                    cutoff = quantile(sjbias, cutoff_quantile);
                    usesj = sjs(sjbias > cutoff(end)); % only take the subjects who are in the highest quantile
                end
                
                alldata = alldata(ismember(alldata.subj_idx, usesj), :);
            end

            % divide RT into quantiles for each subject
            discretizeRTs = @(x) {discretize(x, quantile(x, [0, qntls]))};
            alldata(isnan(alldata.rt), :) = [];

           % ignore if coherence is present but doesn't contain unique values
           if ismember('coherence', alldata.Properties.VariableNames),
           	if length(unique(alldata.coherence(~isnan(alldata.coherence)))) == 1,
           		alldata.coherence = [];
           	end
           end

            % when there were multiple levels of evidence, do these plots
            % separately for each level
            if ~any(ismember('coherence', alldata.Properties.VariableNames))

            	rtbins = splitapply(discretizeRTs, alldata.rt, findgroups(alldata.subj_idx));
            	alldata.rtbins = cat(1, rtbins{:});

                % get RT quantiles for choices that are in line with or against the bias
                [gr, sjidx, rtbins] = findgroups(alldata.subj_idx, alldata.rtbins);
                cpres               = array2table([sjidx, rtbins], 'variablenames', {'subj_idx', 'rtbin'});
                cpres.choice        = splitapply(@nanmean, alldata.biased, gr); % choice proportion
                
                % make into a subjects by rtbin matrix
                mat_tmp = unstack(cpres, 'choice', 'rtbin');
                mat = mat_tmp{:, 2:end}; % remove the last one, only has some weird tail
                
            else
            	disp('splitting by coherence first');

            	[gr, sj, coh] = findgroups(alldata.subj_idx, alldata.coherence);
            	rtbins = splitapply(discretizeRTs, alldata.rt, gr);
            	alldata.rtbins = cat(1, rtbins{:});

                % get RT quantiles for choices that are in line with or against the bias
                [gr, sjidx, rtbins, coh] = findgroups(alldata.subj_idx, alldata.rtbins, alldata.coherence);
                cpres               = array2table([sjidx, rtbins, coh], 'variablenames', {'subj_idx', 'rtbin', 'coh'});
                cpres.choice        = splitapply(@nanmean, alldata.biased, gr); % choice proportion
                
                thesesjs = unique(cpres.subj_idx);
                mat = nan(length(thesesjs), max(cpres.rtbin));
                for sj = 1:length(thesesjs),
                	for r = 1:max(cpres.rtbin);
                		mat(sj, r) = nanmean(cpres.choice(cpres.subj_idx == thesesjs(sj) & cpres.rtbin == r));
                	end
                end
            end
            
            % biased choice proportion
            switch models{m}
            case 'data'
            	disp(size(mat))
                    % ALSO ADD THE REAL DATA WITH SEM/95%CI
                    h = ploterr(qntls, nanmean(mat, 1), [], ...
                    	1.96 * nanstd(mat, [], 1) ./ sqrt(size(mat, 1)), 'k', 'abshhxy', 0);
                    set(h(1), 'color', 'k', 'marker', 'o', ...
                    	'markerfacecolor', 'w', 'markeredgecolor', 'k', 'linewidth', 0.5, 'markersize', 3, ...
                    	'linestyle', '-');
                    set(h(2), 'linewidth', 0.5);
                otherwise
                	plot(qntls, nanmean(mat, 1), 'color', thesecolors{m}, 'linewidth', 1);
                end

            % SAVE         
            if  qidx == 2, % median split, only 2 bins
                % average within each dataset
                avg = splitapply(@nanmean, mat, findgroups(floor(mat_tmp.subj_idx / 1000)));   
                allds.fast(:, m) = avg(:, 1);
                allds.slow(:, m) = avg(:, 2);
            end
        end
        
    axis tight; box off;
    set(gca, 'xtick', qntls);

    axis square;  offsetAxes;
    xlabel('RT (quantiles)');
    set(gca, 'xcolor', 'k', 'ycolor', 'k');
    ylabel(sprintf('P(bias), %s', groups{g}));
    title('Collapsed across datasets');
    set(gca, 'xcolor', 'k', 'ycolor', 'k');

	tightfig;
  	print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/CBFs_collapsed_q%d_sjCutoff%d_%s_models%d.pdf', ...
    qidx, useBiasedSj, groups{g}, whichModels));

%% ========================================== %
% PLOT ACROSS DATASETS - only for median split
% ========================================== %

if qidx ~= 2,
	return;
	disp('not saving');
else
	savefast(sprintf('~/Data/serialHDDM/allds_cbfs_mediansplit.mat'), 'allds');
	disp('saved output');
end
  
end