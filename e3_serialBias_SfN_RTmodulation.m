function e3_serialBias_SfN_RTmodulation

	addpath(genpath('~/code/Tools'));
	warning off; close all; clear;
	global datasets datasetnames colors mypath

	% ========================================== %
	% MODULATION OF SERIAL CHOICE BIAS
	% ========================================== %

	if 0,
	for d = 1:length(datasets),
		close all;
		
		% ========================================== %
		% MODEL WITH ONLY VBIAS
		% ========================================== %
	
		
		% get traces for the model with pupil and rt modulation
		traces = readtable(sprintf('%s/%s/regress_dc_prevresp_prevrt/all_traces.csv', mypath, datasets{d}));

		subplot(4,4,1); hold on;
		[f,xi] = ksdensity(traces.v_prevresp_prevrt);
		h = plot(xi, f, 'color', colors(2, :), 'linewidth', 1);
		% show if these are significant (1-sided?)
		pvalV       = min([mean(traces.v_prevresp_prevrt > 0) mean(traces.v_prevresp_prevrt < 0)]);

		axis tight; axis square;
		title({sprintf('vbias, p = %.3f', pvalV)});
		ylabel('Probability'); vline(0); offsetAxes; xlabel('Parameter (a.u.)');
	
		% ========================================== %
		% MODEL WITH VBIAS AND ZBIAS
		% ========================================== %
	
		subplot(4,4,9); hold on;
						
		% get traces for the model with pupil and rt modulation
		traces = readtable(sprintf('%s/%s/regress_dc_z_prevresp_prevrt/all_traces.csv', mypath, datasets{d}));

		[f,xi] = ksdensity(traces.v_prevresp_prevrt);
		h = plot(xi, f, 'color', colors(2, :), 'linewidth', 1);

		[f,xi] = ksdensity(traces.z_prevresp_prevrt);
		h = plot(xi, f, 'color', colors(1, :), 'linewidth', 1);

		% show if these are significant (1-sided?)
		pvalV       = min([mean(traces.v_prevresp_prevrt > 0) mean(traces.v_prevresp_prevrt < 0)]);
		pvalZ       = min([mean(traces.z_prevresp_prevrt > 0) mean(traces.z_prevresp_prevrt < 0)]);

		axis tight; axis square;
		title({sprintf('vbias, p = %.3f', pvalV); sprintf('zbias, p = %.3f', pvalZ)});
		ylabel('Probability'); vline(0); offsetAxes; xlabel('Parameter (a.u.)');

		% ========================================== %
		% CORRELATION BETWEEN THE TWO
		% ========================================== %

		subplot(4,4,3);

		dat = readtable(sprintf('%s/summary/%s/allindividualresults.csv', ...
		mypath, datasets{d}));
		
		scatter(dat.v_prevresp__regressdcprevrespprevrt, ...
		dat.v_prevrespprevrt__regressdcprevrespprevrt, 10, colors(2, :));
	
		axis square; offsetAxes; xlabel('v_{bias} ~ prevresp'); ylabel('v_{bias} ~ prevresp*prevrt')
		[rho, pval] = corr(dat.v_prevresp__regressdcprevrespprevrt, ...
		dat.v_prevrespprevrt__regressdcprevrespprevrt, 'rows', 'complete');
		title(sprintf('r = %.4f, p = %.4f', rho, pval));
		if pval < 0.05, l = lsline; l.Color = 'k'; 
		end

		subplot(4,4,10);

		dat = readtable(sprintf('%s/summary/%s/allindividualresults.csv', ...
		mypath, datasets{d}));
		
		scatter(dat.z_prevresp__regressdczprevrespprevrt, ...
		dat.z_prevrespprevrt__regressdczprevrespprevrt, 10, colors(1, :));
	
		axis square; offsetAxes; xlabel('z_{bias} ~ prevresp'); ylabel('z_{bias} ~ prevresp*prevrt')
		[rho, pval] = corr(dat.z_prevresp__regressdczprevrespprevrt, ...
		dat.z_prevrespprevrt__regressdczprevrespprevrt, 'rows', 'complete');
		title(sprintf('r = %.4f, p = %.4f', rho, pval));
		if pval < 0.05, l = lsline; l.Color = 'k'; 
		end
		
		subplot(4,4,11);
		
		scatter(dat.v_prevresp__regressdczprevrespprevrt, ...
		dat.v_prevrespprevrt__regressdczprevrespprevrt, 10, colors(2, :));
	
		axis square; offsetAxes; xlabel('v_{bias} ~ prevresp'); ylabel('v_{bias} ~ prevresp*prevrt')
		[rho, pval] = corr(dat.v_prevresp__regressdczprevrespprevrt, ...
		dat.v_prevrespprevrt__regressdczprevrespprevrt, 'rows', 'complete');
		title(sprintf('r = %.4f, p = %.4f', rho, pval));
		if pval < 0.05, l = lsline; l.Color = 'k'; 
		end
    
		% suptitle(datasetnames{d}{1}); 
	    ss = suplabel(cat(2, datasetnames{d}{1}{1}, ' ', datasetnames{d}{1}{2}), 't');
		tightfig;
		
		print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/RT_modulationTraces_d%d.pdf', d));
		
	end
end

	%% ======================================= %%
	% MERGE SUMMARY ACROSS DATASETS
	%% ======================================= %%
	
	close all;
	clear alldat; clear dat;
	for d = 1:length(datasets),
	dat = readtable(sprintf('%s/summary/%s/allindividualresults.csv', ...
	mypath, datasets{d}));
	dat = dat(:, {'subjnr', 'session', 'v_prevresp__regressdcprevrespprevrt', ...
	'v_prevrespprevrt__regressdcprevrespprevrt', 'z_prevresp__regressdczprevrespprevrt', ...
	'z_prevrespprevrt__regressdczprevrespprevrt', 'v_prevresp__regressdczprevrespprevrt', ...
	'v_prevrespprevrt__regressdczprevrespprevrt'});
	
	if d == 1,
		alldat = dat;
	else
		alldat = [alldat; dat];
	end
	end
	dat = alldat; % all observers
	
	subplot(4,4,2);
	scatter(dat.v_prevresp__regressdcprevrespprevrt, ...
	dat.v_prevrespprevrt__regressdcprevrespprevrt, 7, colors(2, :));

	axis square;  axis tight; offsetAxes; 
	% xlabel('v_{bias} ~ prevresp'); 
	ylabel('v_{bias} ~ prevresp*prevrt')
	[rho, pval] = corr(dat.v_prevresp__regressdcprevrespprevrt, ...
	dat.v_prevrespprevrt__regressdcprevrespprevrt, 'rows', 'complete');
	title(sprintf('r = %.4f, p = %.4f', rho, pval));
	if pval < 0.05, l = lsline; l.Color = 'k'; 
	end

	subplot(4,4,5);
	scatter(dat.z_prevresp__regressdczprevrespprevrt, ...
	dat.z_prevrespprevrt__regressdczprevrespprevrt, 7, colors(1, :));

	axis square; axis tight; offsetAxes; xlabel('z_{bias} ~ prevresp'); ylabel('z_{bias} ~ prevresp*prevrt')
	[rho1, pval] = corr(dat.z_prevresp__regressdczprevrespprevrt, ...
	dat.z_prevrespprevrt__regressdczprevrespprevrt, 'rows', 'complete');
	title(sprintf('r = %.4f, p = %.4f', rho1, pval));
	if pval < 0.05, l = lsline; l.Color = 'k'; 
	end
	
	subplot(4,4,6);
	scatter(dat.v_prevresp__regressdczprevrespprevrt, ...
	dat.v_prevrespprevrt__regressdczprevrespprevrt, 7, colors(2, :));

	axis square; axis tight; offsetAxes; xlabel('v_{bias} ~ prevresp'); ylabel('v_{bias} ~ prevresp*prevrt')
	[rho2, pval] = corr(dat.v_prevresp__regressdczprevrespprevrt, ...
	dat.v_prevrespprevrt__regressdczprevrespprevrt, 'rows', 'complete');
	title(sprintf('r = %.4f, p = %.4f', rho2, pval));
	if pval < 0.05, l = lsline; l.Color = 'k'; 
	end
	
	%% steigers test between the two
	[ridiff,~,p] = ridiffci(rho1,rho2,height(dat), height(dat), 0.05);
    ss = suplabel(sprintf('Dr = %.4f, p = %.4f', ridiff, p), 'x');
	% set(ss, 'fontweight', 'normal');
	
	% suptitle(datasetnames{d}{1}); 
    ss = suplabel('All datasets merged', 't');
	% tightfig;
	print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/RT_modulationTraces_merged.pdf'));
	
	
end
