function b2_HDDM_readIntoMatlab(ds, ts, ms)

if ~exist('ds', 'var'), ds =  3:7; end
if ~exist('ts', 'var'), ts = 2;   end
if ~exist('ms', 'var'), ms = 1:9; end

addpath(genpath('~/code/Tools'));
warning off;

usr = getenv('USER');
switch usr
    case 'anne' % local
        datasets = {'RT_RDK', 'projects/0/neurodec/Data/MEG-PL'};
    case 'aeurai' % lisa/cartesius
        datasets = {'RT_RDK', 'MEG-PL', 'MEG-PL-S1', 'MEG-PL-S2', ...
          'Anke-alternating', 'Anke-neutral', 'Anke-repetitive'};
end

for d = ds,
    usepath = sprintf('~/Data/%s/HDDM', datasets{d});
    mdls = {'stimcoding_prevresp_dc', 'stimcoding_prevresp_z', 'stimcoding_prevresp_dc_z', ...
      'regress_dc_prevresp', 'regress_dc_prevresp_prevpupil_prevrt', ...
        'regress_z_prevresp', 'regress_z_prevresp_prevpupil_prevrt', ...
        'regress_dc_z_prevresp', 'regress_dc_z_prevresp_prevpupil_prevrt',};
    traces = {'group', 'all'};

    if ~exist(sprintf('%s/summary', usepath), 'dir'),
      cd(usepath); mkdir('summary');
    end

    switch datasets{d}
        case 'RT_RDK',
            subjects = [3:15 17:25];
        case {'MEG-PL', 'MEG-PL-S1', 'MEG-PL-S2'}
            subjects = 2:65;
        case {'Anke-repetitive', 'Anke-neutral', 'Anke-alternating'}
            subjects = [1:7 9 11:16 18:21 23 24 26 27];
    end

    for m = ms,
        for t = ts,

            chains = 0:29;
            clear dat; clear alldat clear dic;

            % ============================================ %
            % PARSE DIC VALUES FOR THE DIFFERENT CHAINS
            % ============================================ %

            dic = nan(size(chains));
            for c = 1:length(chains),
                fileid = sprintf('%s/%s/DIC-md%d.txt', ...
                    usepath, mdls{m}, chains(c));
                if exist(fileid, 'file'),
                    txtfile = fopen(fileid);
                    txtread = textscan(txtfile, 'Model %d: %d');
                    dic(c)  = txtread{2};
                    fclose(txtfile);
                end
            end

            % ============================================ %
            % CONCATENATE DIFFERENT CHAINS
            % ============================================ %

            for c = 1:length(chains),
                if exist(sprintf('%s/%s/%s_traces-md%d.csv', ...
                        usepath, mdls{m}, traces{t}, chains(c)), 'file'),
                    fprintf('%s/%s/%s_traces-md%d.csv \n', ...
                        usepath, mdls{m}, traces{t}, chains(c));
                    dat{c} = readtable(sprintf('%s/%s/%s_traces-md%d.csv', ...
                        usepath, mdls{m}, traces{t}, chains(c)));
                end
            end

            % skip if these models are not done at all
            if ~exist('dat', 'var'),
              continue;
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
            writetable(alldat, sprintf('%s/summary/%s_%s_traces_concat.csv', ...
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

            tic;
            savefast(sprintf('%s/summary/%s_%s.mat', ...
                usepath, mdls{m}, traces{t}), 'group', 'individuals', 'dic');
            toc;

        end % mdls
    end % traces

    % ============================================ %
    % ONE LARGE TABLE FOR THIS DATASET
    % ============================================ %

    results = array2table(subjects', 'variablenames', {'subjnr'});
    for m = 1:length(mdls),
      try
        disp(mdls{m});
        load(sprintf('%s/summary/%s_%s.mat', usepath, mdls{m}, 'all'));
        disp(mdls{m});
        flds = fieldnames(individuals);
        for p = 1:length(flds),
            if ~isempty(strfind(flds{p}, 'mean')),
                varname = [flds{p}(1:end-5) '__' mdls{m}];
                results.(varname) = individuals.(flds{p});
            end
        end
      end
    end

    writetable(results, sprintf('%s/summary/individualresults.csv', usepath));

end % datasets
end
