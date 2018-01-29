function b2b_Gsq_readIntoMatlab(datasets)

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
    usepath = sprintf('/nfs/aeurai/HDDM/%s/', datasets{d});
    savepath = sprintf('/nfs/aeurai/HDDM/summary/%s', datasets{d});
    if ~exist(savepath, 'dir'),
        cp = pwd; cd('/nfs/aeurai/HDDM/summary');
        mkdir(datasets{d}); cd(cp);
    end
    
    disp(usepath);
    mdls = {'stimcoding_nohist', ...
        'stimcoding_dc_prevresp', ...
        'stimcoding_z_prevresp', ...
        'stimcoding_dc_z_prevresp', ...
        'stimcoding_dc_prevcorrect', ...
        'stimcoding_z_prevcorrect', ...
        'stimcoding_dc_z_prevcorrect', ...
        'regress_dc_z_prevresp_prevstim', ...
        'regress_dc_z_prevresp_prevstim_prevrt_prevpupil', ...
        'regress_dc_z_prevresp_prevstim_prevrt', ...
        'regress_dc_z_prev2resp_prev2stim', ...
        'regress_dc_z_prev3resp_prev3stim', ...
        'stimcoding_nohist_onlyz', 'stimcoding_nohist_onlydc', ...
        'stimcoding_dc_z_prevresp_sessions', ...
        'stimcoding_sz_nohist', ...
        'stimcoding_sz_dc_prevresp', ...
        'stimcoding_sz_z_prevresp', ...
        'stimcoding_sz_dc_z_prevresp', ...
        'stimcoding_dz_z_prevresp_pharma'};
    
    switch datasets{d}
        case 'RT_RDK'
            subjects = [3:15 17:25];
        case {'MEG', 'MEG_MEGsessions'}
            subjects = 2:65;
        case {'Anke_2afc_serial', 'Anke_2afc_neutral', 'Anke_2afc_repetitive', 'Anke_2afc_alternating'},
            subjects = [1:7 9 11:16 18:21 23 24 26 27];
        case 'NatComm'
            subjects = 1:27;
        case {'JW_yesno', 'JW_yesno_2500ms'},
            subjects = [0:23] + 1; % added 1 for matlab indexing
        case 'Anke_merged'
            subjects = [1	2	3	4	5	6	7	9	11	12	13	14	15	16	18	19	20	21	23	24	26	27	31	32	33	34	35	37	38	39	40	41	42	43	44	45];
        case 'Anke_MEG';
            subjects =   [1     2     3     4     5     7     8     9    10    11    12    13    14    15];
        case 'Bharath_fMRI'
            subjects = [4	5	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25];
        case 'Murphy'
            subjects =  [109 110	112	113	114	116	117	119	120	121	124	125	126	127	128	129	131	132	134	135	136	137	138	139	140	142];
    end
    
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
        gsq.bic = BIC.bic;
        
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
    
end % datasets
end
