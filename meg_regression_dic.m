function meg_regression_dic

% Code to fit the history-dependent drift diffusion models described in
% Urai AE, Gee JW de, Donner TH (2018) Choice history biases subsequent evidence accumulation. bioRxiv:251595
%
% MIT License
% Copyright (c) Anne Urai, 2018
% anne.urai@gmail.com


addpath(genpath('~/code/Tools'));
warning off; close all; clear;
global mypath

datasets        = {'MEG_MEGdata'};
datasetnames    = {'MEG trials'};

d = 1;
mdls = {'regress_nohist', ...
    'regress_dc_z_prevresp', ...
    'regress_dc_z_prevresp_motorslope', ...
    'regress_dc_z_motorslope', ...
    'regress_dc_z_prevresp_motorstart', ...
    'regress_dc_z_motorstart', ...
    'regress_dc_z_prevresp_visualgamma', ...
    'regress_dc_z_visualgamma'};

dicvals = [];
for m = 1:length(mdls),
    load(sprintf('%s/summary/%s/%s_all.mat', ...
        mypath, datasets{d}, mdls{m}));
    dicvals = [dicvals mean(dic.chains)];
end

dicvals = dicvals - dicvals(1);
dicvals = dicvals(2:end);
colormap(viridis);
subplot(331);
bar(dicvals);
ylabel({'\Delta DIC from model'; 'without history'}, 'interpreter', 'tex');
box off;
mdlnames = regexprep(mdls(2:end), 'regress_dc_z_', '');
mdlnames = regexprep(mdlnames, '_', ' + ');
set(gca, 'xtick', 1:length(dicvals), ...
    'xticklabel', mdlnames, ...
    'xticklabelrotation',-30);
title('Modulation of dc and z');
tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/MEG_neuralDIC.pdf'));


end
