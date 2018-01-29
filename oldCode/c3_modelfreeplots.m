% from the csv table, make an overview of repetition behaviour

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

set(groot, 'defaultaxesfontsize', 10, 'defaultaxestitlefontsizemultiplier', 1, ...
    'defaultaxestitlefontweight', 'bold', ...
    'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters');

for d = 2:length(datasets),
    
    % load data
    csvfile = dir(sprintf('~/Data/%s/HDDM/*.csv', datasets{d}));
    alldata = readtable(sprintf('~/Data/%s/HDDM/%s', datasets{d}, csvfile.name));
  
    % ============================================ %
    % FIND EXAMPLE PEOPLE
    % ============================================ %
    
    results = readtable(sprintf('~/Data/%s/HDDM/summary/allindividualresults_recoded.csv', datasets{d}));
    

    [val, repeater] = max(results.dc_seq_regress);
    repeater        = results.subjnr(repeater); % get their number
    subjects = unique(results.subjnr)';
    for sj = subjects,
        clf;
        data            = alldata(alldata.subj_idx == sj, :);
        % first, after choice == -1
        trls = find(data.prevresp == -1);
        posterior_predictive(data.stimulus(trls), data.response(trls), data.rt(trls) - 0.5);
        trls = find(data.prevresp == 1);
        posterior_predictive(data.stimulus(trls), data.response(trls), data.rt(trls) - 0.5, 3);
        suplabel('Previous choice 1                Previous choice 0', 'y');
        suplabel(sprintf('P%02d, dc seq = %.3f', sj, results.dc_seq_regress(results.subjnr == sj & results.session == 0)), 't');
        print(gcf, '-dpdf', sprintf('~/Data/%s/HDDM/summary/P%02d_posteriorpredictive.pdf', datasets{d}, sj));
    end
end
