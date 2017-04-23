function e6_serialBias_SfN_modelFree_fromPython

  addpath(genpath('~/code/Tools'));
  warning off; close all; clear;
  set(groot, 'defaultaxesfontsize', 6, 'defaultaxestitlefontsizemultiplier', 1, ...
  'defaultaxestitlefontweight', 'normal', ...
  'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters');

  usr = getenv('USER');
  switch usr
  case 'anne' % local
    datasets = {'RT_RDK', 'projects/0/neurodec/Data/MEG-PL', 'Anke_2afc_sequential', 'NatComm'};
  case 'aeurai' % lisa/cartesius
    datasets = {'RT_RDK', 'MEG', 'Anke_neutral', 'NatComm'};
  end

  params = {'dc', 'z', 'dc_z'};

  % ========================================== %
  % MODELFREE MEASURE OF BIAS
  % RT DISTRIBUTIONS
  % ========================================== %

  for d = 1:4;

    clearvars -except d datasets params datasetnames;
    csvfile = dir(sprintf('~/Data/%s/HDDM/*.csv', datasets{d}));
    csvfile = csvfile(arrayfun(@(x) ~strcmp(x.name(1),'.'), csvfile)); % remove hidden files
    alldata = readtable(sprintf('~/Data/%s/HDDM/%s', datasets{d}, csvfile.name));
    subjects = unique(alldata.subj_idx)';
    disp(subjects);

    for sj = 0:length(subjects)-1,
      clf; cnt = 1;
      for s = 1:2,
        for p = 1:3,
          % this will overlay the two prevresp
          subplot(4,4,cnt); cnt = cnt + 1;
          distmodes(sj+1, s, p, :, :) = plotRTdistributions(datasets{d}, params{p}, sj, s);
          if s == 1, title(regexprep(params{p}, '\_', ' \& '));
          else
            xlabel('Response time (s)');
          end
        end
        cnt = cnt + 1;
      end
      suplabel(sprintf('%s, P%02d', regexprep(datasets{d}, '_', ' '), subjects(sj+1)), 't');
      print(gcf, '-depsc', sprintf('~/Data/serialHDDM/fig4_PostPredPython_%s_P%02d.eps', ...
      datasets{d}, subjects(sj+1)));
      fprintf('~/Data/serialHDDM/fig4_PostPredPython_%s_P%02d.eps \n', ...
      datasets{d}, subjects(sj+1));
    end
    
    savefast(sprintf('~/Data/serialHDDM/fig4_PostPredPython_%s_modes.mat', ...
        datasets{d}), 'distmodes');
    
    % plot distributions of the modes
    % sj, stimulus, param, prevresp, choice
  end
end

function distmodes = plotRTdistributions(d,p,sj,s)
  hold on;

  colors      = cbrewer('div', 'PiYG', 7);
  colors1     = colors([1 end], :);
  colors2     = colors([2 end-1], :);
  edges       = linspace(-3,3,100);
  prevresps   = [-1 1];
  stims       = [0 1];
  distFun     = @(x) histcounts(x, edges, 'normalization', 'pdf');

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % get data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  for c = 1:length(prevresps),
    try
      y = readtable(sprintf('~/Data/%s/HDDM/stimcoding_%s_prevresp/preds/ppq_y_(%d, %d)_subj%d.csv', ...
      d, p, prevresps(c), stims(s), sj), 'readvariablenames', 0);
    catch % for data with coherence levels, take only 10%
      file = dir(sprintf('~/Data/%s/HDDM/stimcoding_%s_prevresp/preds/ppq_y_(0.1*, %.1f, %.1f)_subj%d.csv', ...
      d, p, prevresps(c), stims(s), sj));
      y = readtable(sprintf('~/Data/%s/HDDM/stimcoding_%s_prevresp/preds/%s', d, p, file.name), 'readvariablenames', 0);
    end
    plot(edges, y.Var1, 'color', colors1(c, :), 'linewidth', 1.5);
    [~, maxidx] = max(y.Var1(edges < 0));
    distmodes(c, 1) = edges(maxidx);
    [~, maxidx] = max(y.Var1(edges > 0));
    distmodes(c, 2) = edges(maxidx + length(find(edges < 0)));

  end

  % do separate fits for choice == 1 and choice == -1 (error and correct,
  % depending on the stimulus)
  for c = 1:length(prevresps),
    try
      rt = readtable(sprintf('~/Data/%s/HDDM/stimcoding_%s_prevresp/preds/ppq_rt_(%d, %d)_subj%d.csv', ...
      d, p, prevresps(c), stims(s), sj), 'readvariablenames', 0);
    catch
      file = dir(sprintf('~/Data/%s/HDDM/stimcoding_%s_prevresp/preds/ppq_rt_(0.1*, %.1f, %.1f)_subj%d.csv', ...
      d, p, prevresps(c), stims(s), sj));
      rt = readtable(sprintf('~/Data/%s/HDDM/stimcoding_%s_prevresp/preds/%s', d, p, file.name), 'readvariablenames', 0);
    end
    dist = distFun(rt.Var1);
    plot(edges(1:end-1), dist, 'color', colors2(c, :));
  end

  ylabel('Fraction of trials');
  vline(0);
  text(-2.5, mean(get(gca, 'ylim'))*1.5, 'choice -1', 'fontsize', 5);
  text(1.5, mean(get(gca, 'ylim'))*1.5, 'choice 1', 'fontsize', 5);
  xlim([-3 3]);

end
