addpath(genpath('~/code/Tools'));
warning off; close all; clear;
global datasets datasetnames

usr = getenv('USER');
switch usr
case 'anne' % local
  datasets = {'RT_RDK', 'projects/0/neurodec/Data/MEG-PL', 'NatComm', 'Anke_2afc_neutral'};
case 'aeurai' % lisa/cartesius
  datasets = {'RT_RDK', 'MEG', 'NatComm', 'Anke_neutral', 'Anke_repetitive', 'Anke_alternating'};
end
datasetnames = {'RT', '2IFC',  'NatComm', 'Anke neutral', 'Anke repetitive', 'Anke alternating'};

set(groot, 'defaultaxesfontsize', 7, 'defaultaxestitlefontsizemultiplier', 1.1, ...
    'defaultaxestitlefontweight', 'bold', ...
    'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters');

% ========================================== %
% MODULATION OF SERIAL CHOICE BIAS
% ========================================== %

cnt = 1;
for d = 1:length(datasets),

    dat = readtable(sprintf('~/Data/%s/HDDM/summary/allindividualresults.csv', ...
        datasets{d}));

end

print(gcf, '-depsc', sprintf('~/Data/serialHDDM/fig3_RTmodulation.eps'));
