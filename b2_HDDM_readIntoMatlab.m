function b2_HDDM_readIntoMatlab()

  addpath(genpath('~/code/Tools'));
  close all; clc;
  warning off MATLAB:table:ModifiedVarnames % skip this warning

    datasets = {'RT_RDK', 'MEG', 'MEG_MEGsessions', 'JW_yesno', 'Bharath_fMRI', 'NatComm', 'Anke_2afc_sequential', 'Anke_MEG', 'Anke_merged'};

  for d = 1:length(datasets),
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
    'stimcoding_dc_z_prevresp_sessions'};

    switch datasets{d}
    case 'RT_RDK'
      subjects = [3:15 17:25];
    case {'MEG', 'MEG_MEGsessions'}
      subjects = 2:65;
    case {'Anke_2afc_serial', 'Anke_2afc_neutral', 'Anke_2afc_repetitive', 'Anke_2afc_alternating'},
      subjects = [1:7 9 11:16 18:21 23 24 26 27];
    case 'NatComm'
      subjects = 1:27;
    case 'JW_yesno'
      subjects = [0:23] + 1; % added 1 for matlab indexing
    case 'Anke_merged'
      subjects = [1	2	3	4	5	6	7	9	11	12	13	14	15	16	18	19	20	21	23	24	26	27	31	32	33	34	35	37	38	39	40	41	42	43	44	45];
    case 'Anke_MEG';
      subjects =   [1     2     3     4     5     7     8     9    10    11    12    13    14    15];
    case 'Bharath_fMRI'
      subjects = [4	5	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25];
    end

    for m = 1:length(mdls),

      % skip if this model is empty
      stuff = dir(sprintf('%s/%s', usepath, mdls{m}));
      stuff = stuff(arrayfun(@(x) ~strcmp(x.name(1),'.'),stuff)); % remove hidden stuff
      if isempty(stuff),
        disp(['skipping ' mdls{m}]);
        continue;
      end

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
        assert(all(abs(bsxfun(@minus, dic.chains, dic.full))) < 5, ...
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

        case 'NatComm'

          varnames{v} = regexprep(varnames{v}, '1.0', '1');

          % recode coherence levels
          varnames{v} = regexprep(varnames{v}, '\(0.00625\)', '_c0_0625');
          varnames{v} = regexprep(varnames{v}, '\(0.0125\)', '_c1_25');
          varnames{v} = regexprep(varnames{v}, '\(0.025\)', '_c2_5');
          varnames{v} = regexprep(varnames{v}, '\(0.05\)', '_c5');
          varnames{v} = regexprep(varnames{v}, '\(0.1\)', '_c10');
          varnames{v} = regexprep(varnames{v}, '\(0.2\)', '_c20');
          varnames{v} = regexprep(varnames{v}, '\(0.3\)', '_c30');

        case {'Anke_2afc_neutral', 'Anke_2afc_sequential', 'Anke_2afc_repetitive', 'Anke_2afc_alternating'}

          varnames{v} = regexprep(varnames{v}, '1.0', '1');
          % recode coherence levels in Anke's data
          varnames{v} = regexprep(varnames{v}, '\(0.001\)', '_c0');
          varnames{v} = regexprep(varnames{v}, '\(0.05\)', '_c5');
          varnames{v} = regexprep(varnames{v}, '\(0.1\)', '_c10');
          varnames{v} = regexprep(varnames{v}, '\(0.2\)', '_c20');
          varnames{v} = regexprep(varnames{v}, '\(0.4\)', '_c40');
          varnames{v} = regexprep(varnames{v}, '\(0.6\)', '_c60');

          % some crazy rounding errors
          % recode transition probabilities in Anke's data
          varnames{v} = regexprep(varnames{v}, '\:C\(transitionprob\)\(20.0\)', '_alt');
          varnames{v} = regexprep(varnames{v}, '\:C\(transitionprob\)\(50.0\)', '_neu');
          varnames{v} = regexprep(varnames{v}, '\:C\(transitionprob\)\(80.0\)', '_rep');

          varnames{v} = regexprep(varnames{v}, '20.0\)', '_alt');
          varnames{v} = regexprep(varnames{v}, '50.0\)', '_neu');
          varnames{v} = regexprep(varnames{v}, '80.0\)', '_rep');
        end

        assert(isempty(strfind(varnames{v}, 'transitionprob')), 'no correct parsing')
        % recode some stuff
        varnames{v} = regexprep(varnames{v}, '[().]', '_');
        varnames{v} = regexprep(varnames{v}, '_-1(?=_)', '_2');
        varnames{v} = regexprep(varnames{v}, '_$', '');

        % correct z
        if (isempty(strfind(mdls{m}, 'stimcoding')) & ~isempty(strfind(varnames{v}, 'z_Intercept'))) || ...
          (~isempty(strfind(mdls{m}, 'stimcoding')) & ~isempty(strfind(varnames{v}(1), 'z'))),
          % pointestimates{v, :} = 1 ./ (1 + exp(pointestimates{v, :}));
        end

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
    results = array2table(subjects', 'variablenames', {'subjnr'});
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
            results.(varname) = individuals.(flds{p});
          elseif ~isempty(strfind(flds{p}, 'prct')),
            % also save error bounds around each individual datapoint!
            thismdlname = regexprep(mdls{m}, 'prevresp_prevstim', 'prevrespstim');
            thismdlname = regexprep(thismdlname, 'prevrt_prevpupil', 'prevrtpupil');
            thismdlname = regexprep(thismdlname, 'prevrespstim_prevrtpupil', 'prevrespstimrtpupil');
            thismdlname = regexprep(thismdlname, 'sessions', 'sess');
            thismdlname = regexprep(thismdlname, '_', '');

            % separate variable for the percentiles
            varname = [flds{p}(1:end-5) '_prct__' thismdlname];
            try
              results.(varname) = individuals.(flds{p});
            end
          end
        end
      end
    end

    writetable(results, sprintf('%s/individualresults.csv', savepath));

  end % datasets
end
