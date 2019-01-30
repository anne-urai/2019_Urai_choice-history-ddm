function read_into_Matlab_gSquare(datasets)

% Code to fit the history-dependent drift diffusion models described in
% Urai AE, Gee JW de, Donner TH (2018) Choice history biases subsequent evidence accumulation. bioRxiv:251595
%
% MIT License
% Copyright (c) Anne Urai, 2018
% anne.urai@gmail.com

addpath(genpath('~/code/Tools'));
close all; clc;
warning off MATLAB:table:ModifiedVarnames % skip this warning

for d = 1:length(datasets),

    clear allres;
    mypath  = '/nfs/aeurai/HDDM';
    usepath = sprintf('/nfs/aeurai/HDDM/%s/', datasets{d});
    savepath = sprintf('/nfs/aeurai/HDDM/summary/%s', datasets{d});
    if ~exist(savepath, 'dir'),
        cp = pwd; cd('/nfs/aeurai/HDDM/summary');
        mkdir(datasets{d}); cd(cp);
    end

    disp(usepath);
    % GRAB ALL MODELS THAT EXIST!
    mdls = dir(usepath);
    mdls = mdls(~ismember({mdls.name},{'.','..', '.DS_Store'}));
    mdls = mdls([mdls.isdir]); % only return directories
    mdls = {mdls(:).name};

    % load the csv file to check the existing subject numbers for later
    csvfile = dir(sprintf('%s/%s/*.csv', mypath, datasets{d}));
    csvfile = csvfile(arrayfun(@(x) ~strcmp(x.name(1),'.'), csvfile)); % remove hidden files
    alldata = readtable(sprintf('%s/%s/%s', mypath, datasets{d}, csvfile.name));
    subjects = unique(alldata.subj_idx);  clear alldata;

    for m = 1:length(mdls),

        % skip if this model is empty
        stuff = dir(sprintf('%s/%s', usepath, mdls{m}));
        stuff = stuff(arrayfun(@(x) ~strcmp(x.name(1),'.'),stuff)); % remove hidden files
        if isempty(stuff),
            disp(['skipping ' mdls{m}]);
            continue;
        end

        if ~exist(sprintf('%s/%s/BIC.csv', usepath, mdls{m}), 'file'),
            continue;
        end

        % ============================================ %
        % GET THE GSQUARE FITS
        % ============================================ %

        gsq = readtable(sprintf('%s/%s/Gsquare.csv', usepath, mdls{m}));

        % ============================================ %
        % APPEND BIC VALUES
        % ============================================ %

        BIC = readtable(sprintf('%s/%s/BIC.csv', usepath, mdls{m}));
        % assert(1==0)
        gsq.bic         = BIC.bic;
        gsq.likelihood  = BIC.likelihood;
        gsq.penalty     = BIC.penalty;

        % ============================================ %
        % rename to indicate the model
        % ============================================ %

        gsq.Var1 = [];
        varnames = gsq.Properties.VariableNames;
        varnames(ismember(varnames, 'subj_idx')) = [];
        for f = 1:length(varnames),
            gsq.Properties.VariableNames{varnames{f}} = regexprep([regexprep(varnames{f}, '__1', '_2'), ...
                '__' regexprep(mdls{m}, '_', '')], '___', '__');
        end
        gsq.Properties.VariableNames{'subj_idx'} = 'subjnr';

        if ~exist('allres', 'var')
            allres = gsq;
        else
            allres = join(allres, gsq);
        end
    end

    % remove duplicate columns
    if exist('allres', 'var'),
        writetable(allres, sprintf('%s/individualresults_Gsq.csv', savepath));
    end

    allres.bic__stimcodingnohist

end % datasets
end
