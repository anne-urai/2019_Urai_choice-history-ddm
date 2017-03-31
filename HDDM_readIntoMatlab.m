clear all; close all; warning off;
datasets = {'RT_RDK', 'projects/0/neurodec/Data/MEG-PL'};

for d = 1:length(datasets),
    usepath = sprintf('~/Data/%s/HDDM', datasets{d});
    mdls = {'stimcoding', 'stimcoding_prevresp_dc', 'stimcoding_prevresp_z', ...
        'regress_dc', 'regress_dc_prevresp', 'regress_dc_prevresp_prevpupil_prevrt', ...
        'regress_z_prevresp', 'regress_z_prevresp_prevpupil_prevrt'};
    
    for m = 1:length(mdls),
        chains = 0:59;
        clear dat; clear alldat;
        
        % CONCATENATE DIFFERENT CHAINS
        for c = 1:length(chains),
            dat{c} = readtable(sprintf('%s/%s/group_traces-md%d.csv', ...
                usepath, mdls{m}, chains(c)));
        end
        alldat = cat(1, dat{:});
        
        % rename some variables to make more sense
        for v = 1:length(alldat.Properties.VariableNames),
            
            % make __ into -
            alldat.Properties.VariableNames{v} = ...
                regexprep(alldat.Properties.VariableNames{v}, ...
                '__1_0_', '_0');
            alldat.Properties.VariableNames{v} = ...
                regexprep(alldat.Properties.VariableNames{v}, ...
                '_1_0_', '_1');
            
            % if there were several sessions, remove the underscores around
            alldat.Properties.VariableNames{v} = ...
                regexprep(alldat.Properties.VariableNames{v}, ...
                '_(?!\w)', '');
            
            % forward or backward logistic transform?
            back2prob = @(x) exp(x) ./ (1 + exp(x));
            if strfind(alldat.Properties.VariableNames{v}, 'z_trans'), 
                alldat{:, v} = back2prob(alldat{:, v});
                alldat.Properties.VariableNames{v} = ...
                    regexprep(alldat.Properties.VariableNames{v}, ...
                    'z_trans', 'z');
            end
        end
        
        tic;
        writetable(alldat, sprintf('%s/%s/group_traces_all.csv', ...
            usepath, mdls{m}));
        toc;
        fprintf('%s/%s/group_traces_all.csv \n',  usepath, mdls{m});
    end
end
