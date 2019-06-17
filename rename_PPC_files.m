function b4_renamePPCfiles(datasets)

% Code to fit the history-dependent drift diffusion models as described in
% Urai AE, de Gee JW, Tsetsos K, Donner TH (2019) Choice history biases subsequent evidence accumulation. eLife, in press.
%
% MIT License
% Copyright (c) Anne Urai, 2019
% anne.urai@gmail.com

close all; clc;
mypath  = '/nfs/aeurai/HDDM';
mdls = {'stimcoding_nohist', ...
'stimcoding_dc_prevresp', ...
'stimcoding_z_prevresp', ...
'stimcoding_dc_z_prevresp'};

for d = 1:length(datasets),
for m = 1:length(mdls),
	try
	copyfile(sprintf('%s/%s/%s/ppc_data.csv', mypath, datasets{d}, mdls{m}), ...
	sprintf('%s/summary/%s/%s_ppc_data.csv', mypath, datasets{d}, mdls{m}));
	fprintf('%s/%s/%s/ppc_data.csv \n', mypath, datasets{d}, mdls{m})
end
end
end
