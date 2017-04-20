addpath(genpath('~/code/Tools'));
warning off; close all; clear;
datasets = {'RT_RDK', 'projects/0/neurodec/Data/MEG-PL'};
datasetnames = {'RT', '2IFC', 'Anke'};

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
        dat.v_prevresp__regressdcprevrespstimrtpupilsess, 30, dat.subjnr, '.', ...
        'jitter','on', 'jitterAmount',0.1);
    xlabel('Sessions'); ylabel('v ~ prevresp');
    set(gca, 'xtick', unique(dat.session), 'xlim', [0.5 max(dat.session)]);
    title(datasetnames{d});
    
    subplot(4,4,cnt); cnt = cnt + 1;
    scatter(dat.session, ...
        dat.v_prevrespprevrt__regressdcprevrespstimrtpupilsess, 30, dat.subjnr, '.', ...
        'jitter','on', 'jitterAmount',0.1);
    xlabel('Sessions'); ylabel('v ~ prevresp * prevrt');
    set(gca, 'xtick', unique(dat.session), 'xlim', [0.5 max(dat.session)]);
    title(datasetnames{d});
    
end

print(gcf, '-depsc', sprintf('~/Data/serialHDDM/fig4_stability.eps'));
