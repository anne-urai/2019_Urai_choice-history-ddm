clear all; close all; warning off;
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

  for t = 1:2,
    for m = 1:length(mdls),
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

      % ============================================ %
      % WRITE TO TABLE
      % ============================================ %

      tic;
      writetable(alldat, sprintf('%s/%s/%s_traces_concat.csv', ...
      usepath, mdls{m}, traces{t}));
      toc;
      fprintf('%s/%s/%s_traces_concat.csv \n',  usepath, mdls{m}, traces{t});

      % ============================================ %
      % EXTRACT POINT-ESTIMATES AND QUARTILES
      % ============================================ %

      vars = alldat.Properties.VariableNames;

      for v = 1:length(alldat.Properties.Variablenames),
      end


    end % mdls
  end % traces
end % datasets
