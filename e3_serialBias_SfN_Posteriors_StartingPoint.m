function e3_serialBias_SfN_Posteriors_StartingPoint

	addpath(genpath('~/code/Tools'));
	warning off; close all;
	global mypath datasets datasetnames

	% ========================================== %
	% MODULATION OF SERIAL CHOICE BIAS
	% use the STIMCODING POSTERIORS!
	% ========================================== %

	for d = 1:length(datasets),
		traces = readtable(sprintf('%s/%s/stimcoding_dc_z_prevresp/group_traces.csv', mypath, datasets{d}));
		zbias = invlogit(traces.z_trans_1_) - invlogit(traces.z_trans__1_);
		pval = mean(zbias > 0);
		fprintf('%s %s, mean zbias %.4f, range %.4f to %.4f, p-value %.4f \n', datasetnames{d}{1}{1}, datasetnames{d}{1}{2}, ...
		 mean(zbias), min(zbias), max(zbias), pval);
	end
end
