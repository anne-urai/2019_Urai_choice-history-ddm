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

set(groot, 'defaultaxesfontsize', 4, 'defaultaxestitlefontsizemultiplier', 1, ...
    'defaultaxestitlefontweight', 'bold', ...
    'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters');

for d = 1; %:length(datasets),
    results = readtable(sprintf('~/Data/%s/HDDM/summary/allindividualresults.csv', datasets{d}));
    
    % compute repetition parameters from separate HDDM models
    results.dc_rep__stimcoding_prevresp_dc = ...
        results.dc_1__stimcoding_prevresp_dc - results.dc_2__stimcoding_prevresp_dc;
    
    results.z_rep__stimcoding_prevresp_z = ...
        results.z_2__stimcoding_prevresp_z - results.z_1__stimcoding_prevresp_z;
    
    results.dc_rep__stimcoding_prevresp_dc_z = ...
        results.dc_1__stimcoding_prevresp_dc_z - results.dc_2__stimcoding_prevresp_dc_z;
    
    results.z_rep__stimcoding_prevresp_dc_z = ...
        results.z_2__stimcoding_prevresp_dc_z - results.z_1__stimcoding_prevresp_dc_z;
    
    params = {'repetitioncorr', ...
        'dc_rep__stimcoding_prevresp_dc', 'dc_rep__stimcoding_prevresp_dc_z',...
        'v_prevresp__regress_dc_prevresp','v_prevresp__regress_dc_z_prevresp', ...
        'z_rep__stimcoding_prevresp_z', 'z_rep__stimcoding_prevresp_dc_z', ...
        'z_prevresp__regress_z_prevresp', 'z_prevresp__regress_dc_z_prevresp'};
    
    % see how they correlate
    corrplot(results, params);
    print(gcf, '-dpdf', sprintf('~/Data/%s/HDDM/summary/repetition_corrplot.pdf', datasets{d}));
    
end
