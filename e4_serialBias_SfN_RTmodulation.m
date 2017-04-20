addpath(genpath('~/code/Tools'));
warning off; close all; clear;
global datasets datasetnames
datasets = {'RT_RDK', 'projects/0/neurodec/Data/MEG-PL', 'Anke_2afc_neutral'};
datasetnames = {'RT', '2IFC', 'Anke neutral'};

set(groot, 'defaultaxesfontsize', 7, 'defaultaxestitlefontsizemultiplier', 1.1, ...
    'defaultaxestitlefontweight', 'bold', ...
    'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters');

% ========================================== %
% MODULATION OF SERIAL CHOICE BIAS
% ========================================== %

cnt = 1;
for d = 1:2,
    
    dat = readtable(sprintf('~/Data/%s/HDDM/summary/allindividualresults.csv', ...
        datasets{d}));
    subplot(4,4,cnt);
    scatter(dat.v_prevresp__regressdcprevrespstimrtpupil, ...
        dat.v_prevrespprevrt__regressdcprevrespstimrtpupil, 30, ...
        dat.session, '.');
    lsline;
    xlabel('v ~ prevresp');
    ylabel('v ~ prevresp * prevrt');
    title(datasetnames{d}); axis square;
    
    % also pupil
    subplot(4,4,cnt+4); cnt = cnt + 1;
    scatter(dat.v_prevresp__regressdcprevrespstimrtpupil, ...
        dat.v_prevrespprevpupil__regressdcprevrespstimrtpupil, 30, ...
        dat.session, '.');
    lsline;
    xlabel('v ~ prevresp');
    ylabel('v ~ prevresp * prevpupil');
    title(datasetnames{d}); axis square;
    
end

print(gcf, '-depsc', sprintf('~/Data/serialHDDM/fig3_RTmodulation.eps'));
