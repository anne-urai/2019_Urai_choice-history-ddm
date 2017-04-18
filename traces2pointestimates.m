    if getTraces,
            
            % ============================================ %
            % GET CONCATENATED CHAINS
            % ============================================ %
            
            fprintf('%s/%s/all_traces.csv \n', usepath, mdls{m});
            alldat = readtable(sprintf('%s/%s/all_traces.csv', ...
                usepath, mdls{m}), 'readvariablenames', 1, 'readrownames', 1);
            
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
                back2prob = @(x) 1 ./ (1 + exp(-x));
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
            % EXTRACT POINT ESTIMATES AND QUARTILES
            % ============================================ %
            
            vars        = alldat.Properties.VariableNames';
            paramidx    = find(cellfun(@isempty, strfind(vars, 'subj')) & ...
                cellfun(@isempty, strfind(vars, 'std')) & ...
                cellfun(@isempty, strfind(vars, 'Var1')));
            params      = vars(paramidx); % the parameters of this model
            clear group individuals;
            
            for p = 1:length(params),
                
                % first, group node summary
                group.([params{p} '_mean']) = mean(alldat.(params{p}));
                group.([params{p} '_prct']) = prctile(alldat.(params{p}), [2.5 25 50 75 97.5]);
                
                sjnodes = find(~cellfun(@isempty, strfind(vars, 'subj')) & ...
                    ~cellfun(@isempty, strfind(vars, [params{p} '_'])));
                
                % preallocate subject-specific datapoints
                individuals.([params{p} '_mean']) = nan(fliplr(size(subjects)));
                
                for s = sjnodes',
                    
                    % track the subject specific node for this parameter
                    sjnr    = regexp(vars{s}, '\d+(?=_subj)', 'match'); % numbers at the end
                    sjnr    = str2double(sjnr);
                    sjidx   = find(sjnr == subjects);
                    assert(numel(sjidx) == 1, 'did not find this node');
                    
                    individuals.([params{p} '_mean'])(sjidx)    = mean(alldat.(vars{s}));
                    individuals.([params{p} '_prct'])(sjidx, :) = prctile(alldat.(vars{s}), [2.5 25 50 75 97.5]);
                    
                end
            end
    else
            
        
                
        pointestimates_alldat = varfun(@mean, alldat);
        for v = 1:length(pointestimates_alldat.Properties.VariableNames),
            pointestimates_alldat.Properties.VariableNames{v} = ...
                regexprep(pointestimates_alldat.Properties.VariableNames{v}, ...
                'mean_', '');
        end
        close; hold on;
        for v = 1:length(pointestimates_alldat.Properties.VariableNames),
            varname = pointestimates_alldat.Properties.VariableNames{v};
            plot(pointestimates_alldat.(varname), pointestimates.(varname), 'o');
            % what if they are not the same?
            if abs(pointestimates_alldat.(varname) - pointestimates.(varname)) > 0.001,
                disp(varname);
            end
        end
        waitforbuttonpress;