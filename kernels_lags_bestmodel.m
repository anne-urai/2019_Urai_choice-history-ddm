function kernels_lags_bestmodel

global mypath datasets datasetnames
addpath(genpath('~/code/Tools'));
warning off; close all;

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
    'regress_dcz_lag6', ...
    'regress_z_lag7', ...
    'regress_dc_lag7', ...
    'regress_dcz_lag7'};

alldata.z_correct = nan(length(datasets), 7);
alldata.z_error = nan(length(datasets), 7);
alldata.v_correct = nan(length(datasets), 7);
alldata.v_error = nan(length(datasets), 7);

mat_z.r = nan(length(datasets), 7);
mat_z.pval = nan(length(datasets), 7);
mat_dc.r = nan(length(datasets), 7);
mat_dc.pval = nan(length(datasets), 7);

for d = 1:length(datasets),

	% ============================= %
	% 1. DETERMINE THE BEST MODEL
	% ============================= %

	mdldic = nan(1, length(mdls));
	for m = 1:length(mdls),
	    
	    if ~exist(sprintf('%s/summary/%s/%s_all.mat', ...
	            mypath, datasets{d}, mdls{m}), 'file'),
	        disp('cant find this model')
	        continue;
	    end
	    
	    load(sprintf('%s/summary/%s/%s_all.mat', ...
	        mypath, datasets{d}, mdls{m}));
	    
	    if (isnan(dic.full) || isempty(dic.full)) && ~all(isnan(dic.chains)),
	        dic.full = nanmean(dic.chains);
	    end
	    mdldic(m) = dic.full;
	end

	% everything relative to the full model
	mdldic = bsxfun(@minus, mdldic, mdldic(1));
	mdldic = mdldic(2:end);
	[~, bestMdl] = min(mdldic);
	bestmodelname = regexprep(mdls{bestMdl+1}, '_', '');

	% ========================================================== %
	% 2. FOR THIS MODEL, RECODE INTO CORRECT AND ERROR
	% ========================================================== %

    dat = readtable(sprintf('%s/summary/%s/allindividualresults.csv', mypath, datasets{d}));

	for l = 1:7,
			if l == 1,
				lname = [];
			else
				lname = l;
			end

		try

		alldata.z_correct(d,l) = ...
		        nanmean(dat.(['z_prev' num2str(lname) 'resp__' bestmodelname]) + ...
                dat.(['z_prev' num2str(lname) 'stim__' bestmodelname]));
            end
    end
	for l = 1:7,
					if l == 1,
				lname = [];
			else
				lname = l;
			end

		try
		alldata.z_error(d,l) = ...
		        nanmean(dat.(['z_prev' num2str(lname) 'resp__' bestmodelname]) - ...
                dat.(['z_prev' num2str(lname) 'stim__' bestmodelname]));
            end
    end

	for l = 1:7,
					if l == 1,
				lname = [];
			else
				lname = l;
			end

		try
		alldata.v_correct(d,l) = ...
		         nanmean(dat.(['v_prev' num2str(lname) 'resp__' bestmodelname]) + ...
                dat.(['v_prev' num2str(lname) 'stim__' bestmodelname]));
            end
    end
	for l = 1:7,
			if l == 1,
				lname = [];
			else
				lname = l;
			end

		try
		alldata.v_error(d,l) = ...
		        nanmean(dat.(['v_prev' num2str(lname) 'resp__' bestmodelname]) - ...
                dat.(['v_prev' num2str(lname) 'stim__' bestmodelname]));
            end
    end

    % ALSO COMPUTE CORRELATIONS
    for l = 1:7,

    				if l == 1,
				lname = [];
			else
				lname = l;
			end

try
          [mat_z.r(d, l), mat_z.ci(d,l,:), mat_z.pval(d,l)] = ...
                spearmans(dat.(['z_prev' num2str(lname) 'resp__' bestmodelname]), ...
                dat.(['repetition' num2str(lname)]));
            end
            try
            [mat_dc.r(d, l), mat_dc.ci(d,l,:), mat_dc.pval(d,l)] = ...
                spearmans(dat.(['v_prev' num2str(lname) 'resp__' bestmodelname]), ...
                dat.(['repetition' num2str(lname)]));
            end
 	end

end

    % ========================================================== %
	% 3. PLOT THE VARIABLES THAT ARE PRESENT FOR THIS BEST MODEL
	% ========================================================== %

colors = cbrewer('qual', 'Set2', length(datasets));
models = {'z_correct', 'z_error', 'v_correct', 'v_error'};

% CREATE FIGURE
for pltidx = 1:length(models),

	close all;
	sp1 = subplot(4,4,1); hold on;
	plot([1 7], [0 0], 'k', 'linewidth', 0.5);

	for d = 1:length(datasets),
		plot(1:7, alldata.(models{pltidx})(d, :), 'color', colors(d, :), 'linewidth', 1);
	end
	plot(1:7, nanmean(alldata.(models{pltidx})), 'k', 'linewidth', 1);
	[h, pval] = ttest(alldata.(models{pltidx}));
	if any(h>0),
	    plot(find(h==1), nanmean(alldata.(models{pltidx})(:, (h==1))), ...
	    	'k.', 'markersize', 10);
	end
	xlabel('Lags (# trials)');
	ylabel(regexprep(regexprep(models{pltidx}, '_', ' ~ previous '), 'v ', 'v_{bias} '));
	set(gca, 'xtick', 1:7);
	axis tight; offsetAxes;

	tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/regressionkernels_correcterror_%d.pdf', pltidx));


end


    % ========================================================== %
	% 4. PLOT CORRELATION KERNELS
	% ========================================================== %


% Z CORRELATION KERNELS
close all;
subplot(441); hold on;
plot([1 7], [0 0], 'k', 'linewidth', 0.5);
for d = 1:length(datasets),
    plot(1:7, mat_z.r(d, :), 'color', colors(d, :), 'linewidth', 1);
    h = (mat_z.r(d,:) < 0.05);
    if any(h>0),
        plot(find(h==1), mat_z.r(d, (h==1)), '.', 'markeredgecolor', colors(d, :), ...
        'markerfacecolor', colors(d,:), 'markersize', 5);
    end
end

plot(1:7, nanmean(mat_z.r), 'k', 'linewidth', 1);
[h, pval] = ttest(mat_z.r);
if any(h>0),
    plot(find(h==1), nanmean(mat_z.r(:, (h==1))), ...
    	'.k',  'markersize', 10);
end

ylabel({'Correlation\rho, P(repeat) with' 'z ~ previous response'})
set(gca, 'xtick', 1:7, 'xcolor', 'k', 'ycolor', 'k');
xlabel('Lag (# trials)');
axis square; axis tight; 
ylim([-0.5 1]);
offsetAxes; tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/correlationkernels_z.pdf'));

% DC CORRELATION KERNELS
close all;
subplot(441); hold on;
plot([1 7], [0 0], 'k', 'linewidth', 0.5);
for d = 1:length(datasets),
    plot(1:7, mat_dc.r(d, :), 'color', colors(d, :), 'linewidth', 1);
    h = (mat_dc.pval(d,:) < 0.05);
    if any(h>0),
        plot(find(h==1), mat_dc.r(d, (h==1)), '.', 'markeredgecolor', colors(d, :), ...
        'markerfacecolor', colors(d,:), 'markersize', 5);
    end
end

plot(1:7, nanmean(mat_dc.r), 'k', 'linewidth', 1);
[h, pval] = ttest(mat_dc.r);
if any(h>0),
    plot(find(h==1), nanmean(mat_dc.r(:, (h==1))), ...
    	'.k',  'markersize', 10);
end

ylabel({'Correlation\rho, P(repeat) with' 'v_{bias} ~ previous response'})
set(gca, 'xtick', 1:7, 'xcolor', 'k', 'ycolor', 'k');
xlabel('Lag (# trials)');
axis square; axis tight; 
ylim([-0.5 1]);
offsetAxes; tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/correlationkernels_dc.pdf'));

end
