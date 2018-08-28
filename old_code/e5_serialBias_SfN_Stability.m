
% Code to fit the history-dependent drift diffusion models described in
% Urai AE, Gee JW de, Donner TH (2018) Choice history biases subsequent evidence accumulation. bioRxiv:251595
%
% MIT License
% Copyright (c) Anne Urai, 2018
% anne.urai@gmail.com

addpath(genpath('~/code/Tools'));
warning off; close all; clear;

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
% STABILITY OVER SESSIONS
% ========================================== %

cnt = 1;
for d = 1:2,
    dat = readtable(sprintf('~/Data/%s/HDDM/summary/allindividualresults.csv', ...
        datasets{d}));
    subplot(4,4,cnt); cnt = cnt + 1;
    scatter(dat.session, ...
        dat.v_prevresp__regressdcprevrespstimvasessprevrtpupil, 30, dat.subjnr, '.', ...
        'jitter','on', 'jitterAmount',0.1);
    xlabel('Sessions'); ylabel('v ~ prevresp');
    set(gca, 'xtick', unique(dat.session), 'xlim', [0.5 max(dat.session)]);
    title(datasetnames{d});

    subplot(4,4,cnt); cnt = cnt + 1;
    scatter(dat.session, ...
        dat.v_prevrespprevrt__regressdcprevrespstimvasessprevrtpupil, 30, dat.subjnr, '.', ...
        'jitter','on', 'jitterAmount',0.1);
    xlabel('Sessions'); ylabel('v ~ prevresp * prevrt');
    set(gca, 'xtick', unique(dat.session), 'xlim', [0.5 max(dat.session)]);
    title(datasetnames{d});

end

print(gcf, '-depsc', sprintf('~/Data/serialHDDM/fig4_stability.eps'));
