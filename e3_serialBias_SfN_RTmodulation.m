function e3_serialBias_SfN_RTmodulation

	addpath(genpath('~/code/Tools'));
	warning off; close all; clear;
	global datasets datasetnames colors mypath

	% ========================================== %
	% MODULATION OF SERIAL CHOICE BIAS
	% ========================================== %

	for d = 1:length(datasets),
		close all;
		
		% ========================================== %
		% MODEL WITH ONLY VBIAS
		% ========================================== %
	
		try
			% get traces for the model with pupil and rt modulation
			traces = readtable(sprintf('%s/%s/regress_dc_prevresp_prevrt/all_traces.csv', mypath, datasets{d}));

			subplot(4,4,1); hold on;
			[f,xi] = ksdensity(traces.v_prevresp_prevrt);
			h = plot(xi, f, 'color', colors(2, :), 'linewidth', 1);
			% show if these are significant (1-sided?)
			pvalV       = min([mean(traces.v_prevresp_prevrt > 0) mean(traces.v_prevresp_prevrt < 0)]);
    
			axis tight; axis square;
			title({sprintf('vbias, p = %.3f', pvalV)});
			l.Position(2) = l.Position(2) - 0.15;
			legend boxoff;
			ylabel('Probability');
		end
	
		% ========================================== %
		% MODEL WITH VBIAS AND ZBIAS
		% ========================================== %
	
		subplot(4,4,9); hold on;
			
		try
			
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
			l.Position(2) = l.Position(2) - 0.15;
			legend boxoff;
			ylabel('Probability');
		end

		% ========================================== %
		% CORRELATION BETWEEN THE TWO
		% ========================================== %

		subplot(4,4,3);

		dat = readtable(sprintf('%s/summary/%s/allindividualresults.csv', ...
		mypath, datasets{d}));
		
		try
			scatter(dat.v_prevresp__regressdcprevrespprevrt, ...
			dat.v_prevrespprevrt__regressdcprevrespprevrt, 10, [0.5 0.5 0.5]);
		
			axis square; offsetAxes; xlabel('v_{bias} ~ prevresp'); ylabel('v_{bias} ~ prevresp*prevrt')
			[rho, pval] = corr(dat.v_prevresp__regressdcprevrespprevrt, ...
			dat.v_prevrespprevrt__regressdcprevrespprevrt, 'rows', 'complete');
			title(sprintf('r = %.4f, p = %.4f', rho, pval));
			if pval < 0.05, l = lsline; l.Color = 'k'; 
			end
    
		end
		
		subplot(4,4,10);

		dat = readtable(sprintf('%s/summary/%s/allindividualresults.csv', ...
		mypath, datasets{d}));
		
		try
			scatter(dat.z_prevresp__regressdczprevrespprevrt, ...
			dat.z_prevrespprevrt__regressdczprevrespprevrt, 10, [0.5 0.5 0.5]);
		
			axis square; offsetAxes; xlabel('z_{bias} ~ prevresp'); ylabel('z_{bias} ~ prevresp*prevrt')
			[rho, pval] = corr(dat.z_prevresp__regressdczprevrespprevrt, ...
			dat.z_prevrespprevrt__regressdczprevrespprevrt, 'rows', 'complete');
			title(sprintf('r = %.4f, p = %.4f', rho, pval));
			if pval < 0.05, l = lsline; l.Color = 'k'; 
			end
    
		end
		
		subplot(4,4,11);
		
		try
			scatter(dat.v_prevresp__regressdczprevrespprevrt, ...
			dat.v_prevrespprevrt__regressdczprevrespprevrt, 10, [0.5 0.5 0.5]);
		
			axis square; offsetAxes; xlabel('v_{bias} ~ prevresp'); ylabel('v_{bias} ~ prevresp*prevrt')
			[rho, pval] = corr(dat.v_prevresp__regressdczprevrespprevrt, ...
			dat.v_prevrespprevrt__regressdczprevrespprevrt, 'rows', 'complete');
			title(sprintf('r = %.4f, p = %.4f', rho, pval));
			if pval < 0.05, l = lsline; l.Color = 'k'; 
			end
    
		end
		
	
		% suptitle(datasetnames{d}{1}); 
	    ss = suplabel(cat(2, datasetnames{d}{1}{1}, ' ', datasetnames{d}{1}{2}), 't');
		tightfig;
		
		print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/RT_modulationTraces_d%d.pdf', d));
		
	end
end
