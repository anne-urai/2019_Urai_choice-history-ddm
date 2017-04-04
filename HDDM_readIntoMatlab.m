clear; close all;
warning off;

usr = getenv('USER');
switch usr
    case 'anne' % local
        datasets = {'RT_RDK', 'projects/0/neurodec/Data/MEG-PL'};
    case 'aeurai' % lisa/cartesius
        datasets = {'RT_RDK', 'MEG-PL'};
end

for d = 1:length(datasets),
    usepath = sprintf('~/Data/%s/HDDM', datasets{d});
    mdls = {'stimcoding', 'stimcoding_prevresp_dc', 'stimcoding_prevresp_z', 'stimcoding_prevresp_dc_z', ...
        'regress_dc', 'regress_dc_prevresp', 'regress_dc_prevresp_prevpupil_prevrt', ...
        'regress_z_prevresp', 'regress_z_prevresp_prevpupil_prevrt'};
    traces = {'group', 'all'};
    
    for t = 2,
        for m = 9:length(mdls),
            chains = 0:59;
            clear dat; clear alldat;
            
            % ============================================ %
            % CONCATENATE DIFFERENT CHAINS
            % ============================================ %
            
            for c = 1:length(chains),
                if exist(sprintf('%s/%s/%s_traces-md%d.csv', ...
                        usepath, mdls{m}, traces{t}, chains(c)), 'file'),
                    dat{c} = readtable(sprintf('%s/%s/%s_traces-md%d.csv', ...
                        usepath, mdls{m}, traces{t}, chains(c)));
                end
            end
            alldat = cat(1, dat{:});
            
            % ============================================ %
            % rename some variables to make more sense
            % ============================================ %
            
            for v = 1:length(alldat.Properties.VariableNames),
                
                % remove _0 if preceded by a number
                alldat.Properties.VariableNames{v} = ...
                    regexprep(alldat.Properties.VariableNames{v}, ...
                    '(?<=\d)_0', '');
                
                % recode effect-coded conditions
                alldat.Properties.VariableNames{v} = ...
                    regexprep(alldat.Properties.VariableNames{v}, ...
                    '__1(?=_)', '_2');
                
                % remove extra underscores
                alldat.Properties.VariableNames{v} = ...
                    regexprep(alldat.Properties.VariableNames{v}, ...
                    '__', '_');
                
                % if there were several sessions, remove the underscores around
                alldat.Properties.VariableNames{v} = ...
                    regexprep(alldat.Properties.VariableNames{v}, ...
                    '_(?!\w)', '');
                
                % forward or backward logistic transform? try out!
                back2prob = @(x) 1 ./ (1 + exp(x));
                if strfind(alldat.Properties.VariableNames{v}, '_trans'),
                    alldat{:, v} = back2prob(alldat{:, v});
                    alldat.Properties.VariableNames{v} = ...
                        regexprep(alldat.Properties.VariableNames{v}, ...
                        '_trans', '');
                end
                
                % recode subject nodes
                if strfind(alldat.Properties.VariableNames{v}, '_subj'),
                    % remove this in the middle
                    alldat.Properties.VariableNames{v} = ...
                        regexprep(alldat.Properties.VariableNames{v}, ...
                        '_subj_', '_');
                    alldat.Properties.VariableNames{v} = ...
                        [alldat.Properties.VariableNames{v} '_subj'];
                end
                
            end
            
            % ============================================ %
            % WRITE TO TABLE
            % ============================================ %
            
            tic;
            writetable(alldat, sprintf('%s/%s/%s_traces_concat.csv', ...
                usepath, mdls{m}, traces{t}));
            toc;
            fprintf('%s/%s/%s_traces_concat.csv \n',  usepath, mdls{m}, traces{t});
            
            % ============================================ %
            % EXTRACT POINT ESTIMATES AND QUARTILES
            % ============================================ %
            
            vars        = alldat.Properties.VariableNames';
            paramidx    = find(cellfun(@isempty, strfind(vars, 'subj')) & ...
                cellfun(@isempty, strfind(vars, 'std')) & ...
                cellfun(@isempty, strfind(vars, 'Var1')));
            params      = vars(paramidx); % the parameters of this model
            clear newdat;
            
            for p = 1:length(params),
                
                % first, group node summary
                newdat.([params{p} '_group']).mean = mean(alldat.(params{p}));
                newdat.([params{p} '_group']).prct = prctile(alldat.(params{p}), [2.5 25 50 75 97.5]);
                
                sjnodes = find(~cellfun(@isempty, strfind(vars, 'subj')) & ...
                    ~cellfun(@isempty, strfind(vars, [params{p} '_'])));
                
                switch datasets{d}
                    case 'RT_RDK';
                        subjects = [3:15 17:25];
                end
                
                % preallocate subject-specific datapoints
                newdat.([params{p} '_sj']).mean = nan(fliplr(size(subjects)));
                
                for s = sjnodes',
                    
                    % track the subject specific node for this parameter
                    sjnr    = regexp(vars{s}, '\d+(?=_subj)', 'match'); % numbers at the end
                    sjnr    = str2double(sjnr);
                    sjidx   = find(sjnr == subjects);
                    assert(numel(sjidx) == 1, 'did not find this node');
                    
                    newdat.([params{p} '_sj']).mean(sjidx)    = mean(alldat.(vars{s}));
                    newdat.([params{p} '_sj']).prct(sjidx, :) = prctile(alldat.(vars{s}), [2.5 25 50 75 97.5]);
                    
                end
            end
            
            tic;
            savefast(sprintf('%s/%s/%s_map.mat', ...
                usepath, mdls{m}, traces{t}), 'newdat');
            toc;
            
        end % mdls
    end % traces
end % datasets
