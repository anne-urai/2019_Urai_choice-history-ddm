% get a huge list with values for each participant
% can then work with this dataframe

clear; close all; clc;
addpath(genpath('~/code/Tools'));
warning off;

% ============================================ %
% TWO DIFFERENT DATASETS
% ============================================ %

usr = getenv('USER');
switch usr
    case 'anne' % local
        datasets = {'RT_RDK', 'projects/0/neurodec/Data/MEG-PL'};
    case 'aeurai' % lisa/cartesius
        datasets = {'RT_RDK', 'MEG-PL'};
end

set(groot, 'defaultaxesfontsize', 7, 'defaultaxestitlefontsizemultiplier', 1, ...
    'defaultaxestitlefontweight', 'bold', ...
    'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters');

for d = 1:length(datasets),
    
    % load data
    csvfile = dir(sprintf('~/Data/%s/HDDM/*allsj.csv', datasets{d}));
    alldata = readtable(sprintf('~/Data/%s/HDDM/%s', datasets{d}, csvfile.name));
    
    % rename some things
    try
        alldata.Properties.VariableNames{'trialnr'}     = 'trial';
        alldata.Properties.VariableNames{'blocknr'}     = 'block';
        
        alldata.session = alldata.session + 1; % start at 1
        % session 0 means the average of all sessions
    catch
        alldata.Properties.VariableNames{'subjnr'}   = 'subj_idx';
        alldata.Properties.VariableNames{'stim'}     = 'stimulus';
        alldata.Properties.VariableNames{'resp'}     = 'response';
        alldata.prevrt = circshift(alldata.rt, 1);
    end

    % compute a bunch of basic things from Matlab
    results        = b4_behaviouralMetrics(alldata);
    
end