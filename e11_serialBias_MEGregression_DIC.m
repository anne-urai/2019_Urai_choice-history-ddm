function e11_serialBias_MEGregression_DIC
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
    'regress_dc_z_motorslope', ...
    'regress_dc_z_motorstart', ...
    'regress_dc_z_visualgamma'};

dicvals = [];
for m = 1:length(mdls),
    
    load(sprintf('%s/summary/%s/%s_all.mat', ...
        mypath, datasets{d}, mdls{m}));
    dicvals = [dicvals mean(dic.chains)];
end

dicvals = dicvals - dicvals(1);
dicvals = dicvals(2:end);
bar(dicvals);
ylabel({'\Delta DIC from model'; 'without neural data'}, 'interpreter', 'tex');
box off;
set(gca, 'xtick', 1:length(dicvals), ...
    'xticklabel', {'motor slope', 'motor starting point', 'visual gamma'});

end
