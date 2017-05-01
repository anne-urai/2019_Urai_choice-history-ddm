addpath(genpath('~/code/Tools'));
warning off; close all; clear;
global datasets datasetnames

usr = getenv('USER');
switch usr
case 'anne' % local
  datasets = {'RT_RDK', 'projects/0/neurodec/Data/MEG-PL', 'NatComm', 'Anke_2afc_neutral'};
case 'aeurai' % lisa/cartesius
  datasets = {'RT_RDK', 'MEG', 'NatComm', 'Anke_neutral'};
end
datasetnames = {'RT', '2IFC',  'NatComm', 'Anke neutral', 'Anke repetitive', 'Anke alternating'};

set(groot, 'defaultaxesfontsize', 6, 'defaultaxestitlefontsizemultiplier', 1.1, ...
    'defaultaxestitlefontweight', 'bold', ...
    'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters');

% ========================================== %
% MODULATION OF SERIAL CHOICE BIAS
% ========================================== %

cnt = 1; close all;
for d = 1:length(datasets),

  % get traces for the model with pupil and rt modulation
  traces = readtable(sprintf('~/Data/%s/HDDM/regress_dc_prevresp_prevstim_prevrt_prevpupil/all_traces.csv', datasets{d}));

  % plot the pupil and RT traces
  subplot(4,4,d); hold on;
  histogram(traces.v_prevresp_prevrt, 'displaystyle', 'stairs');
  histogram(traces.v_prevresp_prevpupil, 'displaystyle', 'stairs');

  % show if these are significant (1-sided?)
  pvalRT = mean(traces.v_prevresp_prevrt > 0);
  pvalPupil = mean(traces.v_prevresp_prevpupil > 0);

  axis tight;
  l = legend({sprintf('RT, p = %.3f', pvalRT); sprintf('Pupil, p = %.3f', pvalPupil)}, ...
  'location', 'southeast');
  l.Position(2) = l.Position(2) - 0.15;
  legend boxoff;
  title(datasetnames{d}); xlabel('Prevresp modulation');
  vline(0);

end

print(gcf, '-depsc', sprintf('~/Data/serialHDDM/fig3_RTmodulation_traces.eps'));
