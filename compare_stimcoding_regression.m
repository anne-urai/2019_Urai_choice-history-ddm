function barplots_modelcomparison_regression()

% Code to fit the history-dependent drift diffusion models described in
% Urai AE, Gee JW de, Donner TH (2018) Choice history biases subsequent evidence accumulation. bioRxiv:251595
%
% MIT License
% Copyright (c) Anne Urai, 2018
% anne.urai@gmail.com

addpath(genpath('~/code/Tools'));
warning off; close all;
global datasets datasetnames mypath colors

for d = 1:length(datasets),

    try
    dat = readtable(sprintf('%s/summary/%s/allindividualresults.csv', mypath, datasets{d}));
    dat = dat(dat.session == 0, :);

    whichmdls = 'stimcoding';
       dat.z_corr = ...
        dat.(['z_1_1__' whichmdls 'dczprevcorrect']) - dat.(['z_1_2__' whichmdls 'dczprevcorrect']);
    dat.v_corr = ...
        dat.(['dc_1_1__' whichmdls 'dczprevcorrect']) - dat.(['dc_1_2__' whichmdls 'dczprevcorrect']);
    
        dat.z_err = ...
            dat.(['z_0_1__' whichmdls 'dczprevcorrect']) - dat.(['z_0_2__' whichmdls 'dczprevcorrect']);
        dat.v_err = ...
            dat.(['dc_0_1__' whichmdls 'dczprevcorrect']) - dat.(['dc_0_2__' whichmdls 'dczprevcorrect']);
 

    fullmodelname = 'regressdczlag1';

    dat.z_corrR = dat.(['z_prev' 'resp__' fullmodelname]) + ...
                    dat.(['z_prev' 'stim__' fullmodelname]);
    dat.z_errR = dat.(['z_prev' 'resp__' fullmodelname]) - ...
                    dat.(['z_prev' 'stim__' fullmodelname]);
    dat.v_corrR = dat.(['v_prev' 'resp__' fullmodelname]) + ...
                    dat.(['v_prev' 'stim__' fullmodelname]);
    dat.v_errR = dat.(['v_prev' 'resp__' fullmodelname]) - ...
                    dat.(['v_prev' 'stim__' fullmodelname]);
    close all;
    corrplot(dat, {'z_corr', 'z_err', 'v_corr', 'v_err'}, ...
        {'z_corrR', 'z_errR', 'v_corrR', 'v_errR'});

    suptitle(cat(2, datasetnames{d}{1}, ' ', datasetnames{d}{2}));
    %title(datasetnames);
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/stimcoding_vs_regression_d%d.pdf', d));
    fprintf('~/Data/serialHDDM/stimcoding_vs_regression_d%d.pdf \n', d);
end
end