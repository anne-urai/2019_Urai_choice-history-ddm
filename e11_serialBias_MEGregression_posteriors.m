function e11_serialBias_MEGregression_posteriors
% ========================================== %
% conditional response functions from White & Poldrack
% run on simulated rather than real data
% ========================================== %

addpath(genpath('~/code/Tools'));
warning off; close all; clear;
global mypath

datasets        = {'MEG_MEGdata'};
datasetnames    = {'MEG trials'};
d = 1;

mdls = {'regress_dc_z_prevresp', ...
    'regress_dc_z_prevresp_motorslope', ...
    'regress_dc_z_motorslope', ...
    'regress_dc_z_prevresp_motorstart', ...
    'regress_dc_z_motorstart', ...
    'regress_dc_z_prevresp_visualgamma', ...
    'regress_dc_z_visualgamma'};

for m = 1:length(mdls),
    close all;
    dat = readtable(sprintf('%s/%s/%s/group_traces.csv', mypath, datasets{d}, mdls{m}));
    
    if strfind(mdls{m}, 'motorstart'),
        dat.v_motorstart = dat.v_motorbeta;
        dat.z_motorstart = dat.z_motorbeta;
    end
    
    % which parameters do we need?
    params = regexprep(mdls{m}, 'regress_dc_z_', '');
    params = strsplit(params, '_');
    for p = 1:length(params),
        subplot(4,4,(p-1)*4+1); hold on; plotHist(dat, ['v_' params{p}]);
        subplot(4,4,(p-1)*4+2); hold on; plotHist(dat, ['z_' params{p}]);
    end
    suplabel('Posterior probability', 'y');
    tightfig;
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/MEG_posteriors_%s.pdf', mdls{m}));
end

end

function plotHist(dat, param)
hold on; histogram(dat.(param), 'edgecolor', 'none', 'facecolor', [0.5 0.5 0.5]); vline(0, 'color', 'r');
xlabel({regexprep(param, '_', '~'); sprintf('p = %.4f', min([mean(dat.(param) < 0) mean(dat.(param) > 0)]))});
end