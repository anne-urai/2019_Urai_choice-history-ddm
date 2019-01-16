function b2_HDDM_readIntoMatlab(datasets)

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
    
    usr = getenv('USER');
    switch usr
        case 'anne'
            mypath = '~/Data/HDDM';
        case 'aeurai'
            mypath  = '/nfs/aeurai/HDDM';
    end
    
    usepath = sprintf('%s/%s/', mypath, datasets{d});
    savepath = sprintf('%s/summary/%s', mypath, datasets{d});
    if ~exist(savepath, 'dir'),
        cp = pwd; cd(sprintf('%s/summary', mypath));
        mkdir(datasets{d}); cd(cp);
    end
    
    disp(usepath);
    mdls = {'stimcoding_nohist', ...
        'stimcoding_dc_prevresp', ...
        'stimcoding_z_prevresp', ...
        'stimcoding_dc_z_prevresp', ...
        'stimcoding_dc_z_prevresp_st', ...
        'stimcoding_dc_z_prevresp_pharma', ...
        'stimcoding_dc_z_prevcorrect', ...
        'stimcoding_prevcorrect', ...
        'regress_nohist', ...
        'regress_dc_lag1', ...
        'regress_z_lag1', ...
        'regress_dcz_lag1', ...
        'regress_dc_lag2', ...
        'regress_z_lag2', ...
        'regress_dcz_lag2', ...
        'regress_dc_lag3', ...
        'regress_z_lag3', ...
        'regress_dcz_lag3', ...
        'stimcoding_nohist_stcoh',  ...
        'stimcoding_dc_prevresp_stcoh',  ...
        'stimcoding_z_prevresp_stcoh',  ...
        'stimcoding_dc_z_prevresp_stcoh'};
    
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
        
        chains = 0:29;
        clear dat alldat dic;
        
        % ============================================ %
        % PARSE DIC VALUES FOR THE DIFFERENT CHAINS
        % ============================================ %
        
        dic.chains = nan(size(chains));
        for c = 1:length(chains),
            fileid = sprintf('%s/%s/DIC-md%d.txt', ...
                usepath, mdls{m}, chains(c));
            if exist(fileid, 'file'),
                txtfile = fopen(fileid);
                txtread = textscan(txtfile, 'Model %d: %.3f');
                dic.chains(c)  = txtread{2};
                fclose(txtfile);
            end
        end
        
        % ALSO GET THE ONE FROM THE CONCATENATED MODEL
        fileid = sprintf('%s/%s/DIC-combined.txt', ...
            usepath, mdls{m});
        if exist(fileid, 'file'),
            txtfile = fopen(fileid);
            txtread = textscan(txtfile, 'Combined model: %.3f');
            dic.full = double(txtread{1});
            if dic.full == 0, dic.full = NaN;
            end
            fclose(txtfile);
        else
            dic.full = NaN;
        end
        
        if ~isempty(dic.full),
            % CHECK THAT THE DICS HAVE CONVERGED
            assert(all(abs(bsxfun(@minus, dic.chains, dic.full))) < 1.05, ...
                'chains differ in DIC');
        else
            assert(all(abs(bsxfun(@minus, dic.chains, nanmean(dic.chains)))) < 5, ...
                'chains differ in DIC');
        end
        
        % ============================================ %
        % CHECK R-HAT GELMAN-RUBIN STATISTIC
        % ============================================ %
        
        try
            rhat = readtable(sprintf('%s/%s/gelman_rubin.txt', usepath, mdls{m}));
            notConverged = rhat(find(abs(rhat.Var2 - 1) > 0.001), :);
            % ignore single-subject convergence or std
            notConverged(find(~cellfun(@isempty, strfind(notConverged.Var1, 'subj'))), :) = [];
            notConverged(find(~cellfun(@isempty, strfind(notConverged.Var1, 'std'))), :) = [];
            if ~isempty(notConverged),
                % warning('not all group parameters have converged');
                %  disp(notConverged);
            end
        end
        
        % ============================================ %
        % ALSO GET POINT ESTIMATES FROM RESULTS FILE
        % ============================================ %
        
        if ~exist(sprintf('%s/%s/results-combined.csv', ...
                usepath, mdls{m}), 'file'),
            disp('skipping');
            continue;
        end
        
        % compare with results
        pointestimates = readtable(sprintf('%s/%s/results-combined.csv', ...
            usepath, mdls{m}), 'readrownames', 1);
        
        varnames     = pointestimates.Properties.RowNames;
        varnamesOrig = pointestimates.Properties.RowNames;
        
        for v = 1:length(varnames),
            
            % remove trailing .0 for subj nr
            varnames{v} = regexprep(varnames{v}, '\.0$', '');
            varnames{v} = regexprep(varnames{v}, '\[', '(');
            varnames{v} = regexprep(varnames{v}, '\]', ')');
            
            switch datasets{d},
                
                case {'NatComm', 'NatComm_500ms'}
                    
                    % varnames{v} = regexprep(varnames{v}, '0.1', '0_1');
                    varnames{v} = regexprep(varnames{v}, '1.0', '1');
                    
                    % recode coherence levels
                    if ~isempty(strfind(varnames{v}, 'v')),
                        varnames{v} = regexprep(varnames{v}, '\(0.00625\)', '_c0_0625');
                        varnames{v} = regexprep(varnames{v}, '\(0.0125\)', '_c1_25');
                        varnames{v} = regexprep(varnames{v}, '\(0.025\)', '_c2_5');
                        varnames{v} = regexprep(varnames{v}, '\(0.05\)', '_c5');
                        varnames{v} = regexprep(varnames{v}, '\(0.1\)', '_c10');
                        varnames{v} = regexprep(varnames{v}, '\(0.2\)', '_c20');
                        varnames{v} = regexprep(varnames{v}, '\(0.3\)', '_c30');
                    end
                    
                case {'Anke_2afc_neutral', 'Anke_2afc_sequential', 'Anke_2afc_repetitive', 'Anke_2afc_alternating', ...
                        'Bharath_fMRI', 'Anke_MEG', 'Anke_MEG_neutral'}
                    
                    varnames{v} = regexprep(varnames{v}, '1.0', '1');
                    % recode coherence levels in Anke's data
                    varnames{v} = regexprep(varnames{v}, '\(0.001\)', '_c0');
                    varnames{v} = regexprep(varnames{v}, '\(0.05\)', '_c5');
                    varnames{v} = regexprep(varnames{v}, '\(0.1\)', '_c10');
                    varnames{v} = regexprep(varnames{v}, '\(0.2\)', '_c20');
                    varnames{v} = regexprep(varnames{v}, '\(0.4\)', '_c40');
                    varnames{v} = regexprep(varnames{v}, '\(0.6\)', '_c60');
                    
                    % also for anke's meg data
                    varnames{v} = regexprep(varnames{v}, '\(0.0\)', '_c0');
                    varnames{v} = regexprep(varnames{v}, '\(0.03\)', '_c3');
                    varnames{v} = regexprep(varnames{v}, '\(0.09\)', '_c9');
                    varnames{v} = regexprep(varnames{v}, '\(0.81\)', '_c81');
                    varnames{v} = regexprep(varnames{v}, '\(0.27\)', '_c27');
                    
            end
            
            if strfind(varnames{v}, 'transitionprob'),
                
                if strfind(mdls{m}, 'regress'),
                    
                    % rename variables in regression models
                    varnames{v} = regexprep(varnames{v}, 'C\(transitionprob\)\[20.0\]\:prevresp', 'prevresp_alt');
                    varnames{v} = regexprep(varnames{v}, 'C\(transitionprob\)\[50.0\]\:prevresp', 'prevresp_neu');
                    varnames{v} = regexprep(varnames{v}, 'C\(transitionprob\)\[80.0\]\:prevresp', 'prevresp_rep');
                    
                    varnames{v} = regexprep(varnames{v}, 'C\(transitionprob\)\(20.0\)', '_alt');
                    varnames{v} = regexprep(varnames{v}, 'C\(transitionprob\)\(50.0\)', '_neu');
                    varnames{v} = regexprep(varnames{v}, 'C\(transitionprob\)\(80.0\)', '_rep');
                    
                elseif strfind(mdls{m}, 'stimcoding'),
                    
                    % some crazy rounding errors
                    % recode transition probabilities in Anke's data
                    varnames{v} = regexprep(varnames{v}, '\:C\(transitionprob\)\(20.0\)', '_alt');
                    varnames{v} = regexprep(varnames{v}, '\:C\(transitionprob\)\(50.0\)', '_neu');
                    varnames{v} = regexprep(varnames{v}, '\:C\(transitionprob\)\(80.0\)', '_rep');
                    
                    varnames{v} = regexprep(varnames{v}, '20.0\)', '_alt');
                    varnames{v} = regexprep(varnames{v}, '50.0\)', '_neu');
                    varnames{v} = regexprep(varnames{v}, '80.0\)', '_rep');
                end
            end
            
            assert(isempty(strfind(varnames{v}, 'transitionprob')), 'no correct parsing')
            % recode some stuff
            varnames{v} = regexprep(varnames{v}, '[().]', '_');
            varnames{v} = regexprep(varnames{v}, '_-1(?=_)', '_2');
            varnames{v} = regexprep(varnames{v}, '_$', '');
            
            if strfind(varnames{v}, 'session'),
                varnames{v} = regexprep(varnames{v}, 'C[a-zA-Z_]+(?=\d)', '_s');
                varnames{v} = regexprep(varnames{v}, '\]', '');
            end
            
            if strfind(varnames{v}, '_subj') & strfind(mdls{m}, 'stimcoding'),
                % remove this in the middle
                varnames{v} = regexprep(varnames{v},'_subj_', '_');
                sjnum       = regexp(varnames{v}, '\d+$', 'match');
                varnames{v} = regexprep(varnames{v}, '\d+$', '');
                varnames{v} = [varnames{v} '_subj_' sjnum{1}];
            end
            
            varnames{v} = regexprep(varnames{v}, '__', '_');
            varnames{v} = regexprep(varnames{v}, ':', '');
            
        end
        
        pointestimates.Properties.RowNames = varnames;
        vars                = pointestimates.Properties.RowNames;
        paramidx            = find(cellfun(@isempty, strfind(vars, 'subj')) & ...
            cellfun(@isempty, strfind(vars, 'std')) & ...
            cellfun(@isempty, strfind(vars, 'Var1')));
        params              = vars(paramidx); % the parameters of this model
        clear group individuals;
        
        for p = 1:length(params),
            
            % first, group node summary
            group.([params{p} '_mean']) = pointestimates{params{p}, {'mean'}};
            group.([params{p} '_prct']) = pointestimates{params{p}, {'x2_5q', 'x25q', 'x50q', 'x75q', 'x97_5q'}};
            
            sjnodes = find(~cellfun(@isempty, strfind(vars, 'subj')) & ...
                ~cellfun(@isempty, strfind(vars, [params{p} '_'])));
            
            % preallocate subject-specific datapoints
            individuals.([params{p} '_mean']) = nan(fliplr(size(subjects)));
            
            for s = sjnodes',
                
                % track the subject specific node for this parameter
                sjnr    = regexp(vars{s}, '(?<=subj_)\d+', 'match'); % numbers at the end
                sjnr    = str2double(sjnr);
                sjidx   = find(sjnr == subjects);
                assert(numel(sjidx) == 1, 'did not find this node');
                
                % keep the mean and percentiles
                individuals.([params{p} '_mean'])(sjidx)    = pointestimates{vars(s), {'mean'}};
                individuals.([params{p} '_prct'])(sjidx, :) = pointestimates{vars{s}, {'x2_5q', 'x25q', 'x50q', 'x75q', 'x97_5q'}};
                
            end
        end
        
        %tic;
        savefast(sprintf('%s/%s_all.mat', ...
            savepath, mdls{m}), 'group', 'individuals', 'dic');
        %toc;
        
    end % mdls
    
    % ============================================ %
    % ONE LARGE TABLE FOR THIS DATASET
    % ============================================ %
    
    disp('making table');
    results = array2table(subjects, 'variablenames', {'subjnr'});
    results_group = array2table(nan(1,1), 'variablenames', {'subjnr'});
    for m = 1:length(mdls),
        if exist(sprintf('%s/%s_%s.mat', savepath, mdls{m}, 'all'), 'file'),
            load(sprintf('%s/%s_%s.mat', savepath, mdls{m}, 'all'));
            flds = fieldnames(individuals);
            for p = 1:length(flds),
                if ~isempty(strfind(flds{p}, 'mean')),
                    % avoid that names get too long, remove spaces from model name
                    thismdlname = regexprep(mdls{m}, 'prevresp_prevstim', 'prevrespstim');
                    thismdlname = regexprep(thismdlname, 'prevrt_prevpupil', 'prevrtpupil');
                    thismdlname = regexprep(thismdlname, 'prevrespstim_prevrtpupil', 'prevrespstimrtpupil');
                    thismdlname = regexprep(thismdlname, 'sessions', 'sess');
                    thismdlname = regexprep(thismdlname, '_', '');
                    varname = [flds{p}(1:end-5) '__' thismdlname];
                    results.(varname) = individuals.(flds{p})';
                    
                    % also save the group-level
                    results_group.(varname) = group.(flds{p});
                    
                elseif ~isempty(strfind(flds{p}, 'prct')),
                    % also save error bounds around each individual datapoint!
                    thismdlname = regexprep(mdls{m}, 'prevresp_prevstim', 'prevrespstim');
                    thismdlname = regexprep(thismdlname, 'prevrt_prevpupil', 'prevrtpupil');
                    thismdlname = regexprep(thismdlname, 'prevrespstim_prevrtpupil', 'prevrespstimrtpupil');
                    thismdlname = regexprep(thismdlname, 'sessions', 'sess');
                    thismdlname = regexprep(thismdlname, '_', '');
                    
                    % separate variable for the percentiles
                    varname = [flds{p}(1:end-5) '_prct__' thismdlname];
                end
            end
        end
    end
    
    
    % ================================================= %
    % ADD GROUP-LEVEL SUMMARY STATS AND DIC
    % ================================================= %
    
    results.dic = nan(size(results.subjnr));
    results_group.dic = dic.full;
    
    % append
    results_full = [results; results_group];
    % allstuff = load(sprintf('%s/%s_all.mat', savepath, mdls{m}));
    
    writetable(results_full, sprintf('%s/individualresults.csv', savepath));
    
end % datasets
end
