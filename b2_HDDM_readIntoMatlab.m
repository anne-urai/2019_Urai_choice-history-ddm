function b2_HDDM_readIntoMatlab(ds, ts, ms)

if ~exist('ds', 'var'), ds = 1:3; end
if ~exist('ms', 'var'), ms = 1:17; end

addpath(genpath('~/code/Tools'));
warning verbose;
close all; clc;

getTraces = false; % ca
usr = getenv('USER');
switch usr
    case 'anne' % local
        datasets = {'RT_RDK', 'projects/0/neurodec/Data/MEG-PL', 'Anke_2afc_sequential'};
    case 'aeurai' % lisa/cartesius
        datasets = {'RT_RDK', 'MEG', 'Anke_serial'};
end

for d = ds,
    usepath = sprintf('~/Data/%s/HDDM', datasets{d});
    mdls = {'stimcoding_dc_prevresp_prevstim', ...
        'stimcoding_z_prevresp_prevstim', ...
        'stimcoding_dc_z_prevresp_prevstim', ...
        'regress_dc_prevresp_prevstim', ...
        'regress_z_prevresp_prevstim', ...
        'regress_dc_z_prevresp_prevstim', ...
        'regress_dc_prevresp_prevstim_sessions', ...
        'regress_dc_prevresp_prevstim_prevrt', ...
        'regress_dc_prevresp_prevstim_prevrt_sessions', ...
        'regress_dc_prevresp_prevstim_prevrt_prevpupil', ...
        'regress_dc_prevresp_prevstim_prevrt_prevpupil_sessions', ...
        'stimcoding_dc_prevresp', ...
        'stimcoding_z_prevresp', ...
        'stimcoding_dc_z_prevresp', ...
        'regress_dc_prevresp', ...
        'regress_z_prevresp', ...
        'regress_dc_z_prevresp'};
    
    if ~exist(sprintf('%s/summary', usepath), 'dir'),
        cd(usepath); mkdir('summary');
    end
    
    switch d
        case 1 % RT-RDK
            subjects = [3:15 17:25];
        case 2 % MEG
            subjects = 2:65;
        case 3 % Anke
            subjects = [1:7 9 11:16 18:21 23 24 26 27];
    end
    
    for m = ms,
        
        try
            disp(mdls{m});
            chains = 0:14;
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
                txtread = textscan(txtfile, 'Model %d: %.3f');
                dic.full = double(txtread{2});
                if dic.full == 0, dic.full = NaN; end
                fclose(txtfile);
            end
            
            try
                plot(bsxfun(@minus, dic.chains, dic.full), 'o');
                hline(mean(dic.chains) - dic.full, 'b');
                xlabel('Chains'); ylabel('DIC');
                if abs(mean(dic.chains) - dic.full) < 1, ...
                        warning('DIC values differ between chains');
                end
            end
            % waitforbuttonpress; close;
            
            % ============================================ %
            % ALSO GET POINT ESTIMATES FROM RESULTS FILE
            % ============================================ %
            
            % compare with results
            pointestimates = readtable(sprintf('%s/%s/results-combined.csv', ...
                usepath, mdls{m}), 'readrownames', 1);
            
            varnames = pointestimates.Properties.RowNames;
            for v = 1:length(varnames),
                
                % recode some stuff
                varnames{v} = regexprep(varnames{v}, '[().]', '_');
                varnames{v} = regexprep(varnames{v}, '_-1(?=_)', '_2');
                varnames{v} = regexprep(varnames{v}, '_$', '');
                
                if strfind(varnames{v}, 'session'),
                    varnames{v} = regexprep(varnames{v}, 'C[a-zA-Z_0-9\f\n\r\t\v\[]+(?=\d)', '_s');
                    varnames{v} = regexprep(varnames{v}, '\]', '');
                end
                
                if strfind(varnames{v}, '_subj'),
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
                    
                    individuals.([params{p} '_mean'])(sjidx)    = pointestimates{vars(s), {'mean'}};
                    individuals.([params{p} '_prct'])(sjidx, :) = pointestimates{vars{s}, {'x2_5q', 'x25q', 'x50q', 'x75q', 'x97_5q'}};
                    
                end
            end
            
            tic;
            savefast(sprintf('%s/summary/%s_group.mat', ...
                usepath, mdls{m}), 'group', 'individuals', 'dic');
            toc;
            
        end % mdls
    end
    
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
