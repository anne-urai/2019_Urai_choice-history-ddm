function trialCount
% ========================================== %
% conditional response functions from White & Poldrack
% run on simulated rather than real data
% ========================================== %

addpath(genpath('~/code/Tools'));
warning off; close all; clear;
global datasets datasetnames mypath colors

for d = 1:length(datasets),
    
    filename = dir(sprintf('%s/%s/*.csv', mypath, datasets{d}));
    alldata  = readtable(sprintf('%s/%s/%s', mypath, datasets{d}, filename.name));
    
    nrsj(d) = numel(unique(alldata.subj_idx));
    nrtrials(d) = size(alldata, 1);
end