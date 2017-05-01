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
for d = 1:4; %length(datasets),

    dat = readtable(sprintf('~/Data/%s/HDDM/summary/allindividualresults.csv', ...
        datasets{d}));

    % also load the full file for this model
    load(sprintf('~/Data/%s/HDDM/summary/regress_dc_prevresp_prevstim_prevrt_prevpupil_all.mat', ...
        datasets{d}));

    subplot(4,4,cnt); hold on;
    scatter(dat.v_prevresp__regressdcprevrespstimrtpupil, ...
        dat.v_prevrespprevrt__regressdcprevrespstimrtpupil, 30, ...
        dat.session, '.');
    [rho, pval] = corr(dat.v_prevresp__regressdcprevrespstimrtpupil, ...
        dat.v_prevrespprevrt__regressdcprevrespstimrtpupil, ...
        'type', 'spearman', 'rows', 'complete');
    if pval < 0.05, lsline; end;
    xlabel('v ~ prevresp');
    ylabel('v ~ prevresp * prevrt');
    title(datasetnames{d}); axis square;
  %  axis tight;

if 0,
    % add the quantiles on the axes
    sz = [5 10 15 10 5];
    xlm = min(get(gca, 'xlim'));
    ylm = min(get(gca, 'ylim'));

    scatter(ones(1,5) * xlm, group.v_prevrespprevrt_prct, ...
        sz, 'k', '>');
    scatter(group.v_prevresp_prct, ones(1,5) * ylm, ...
        sz, 'k', '^');
    axis tight;
end

    % also pupil
    subplot(4,4,cnt+4); cnt = cnt + 1; hold on;
    scatter(dat.v_prevresp__regressdcprevrespstimrtpupil, ...
        dat.v_prevrespprevpupil__regressdcprevrespstimrtpupil, 30, ...
        dat.session, '.');
    [rho, pval] = corr(dat.v_prevresp__regressdcprevrespstimrtpupil, ...
        dat.v_prevrespprevpupil__regressdcprevrespstimrtpupil, ...
        'type', 'spearman', 'rows', 'complete');
    if pval < 0.05, lsline; end;
    xlabel('v ~ prevresp');
    ylabel('v ~ prevresp * prevpupil');
    title(datasetnames{d}); axis square;
    %axis tight;

if 0,
    % add the quantiles on the axes
    xlm = min(get(gca, 'xlim'));
    ylm = min(get(gca, 'ylim'));
    scatter(ones(1,5) * xlm, group.v_prevrespprevpupil_prct, ...
        sz, 'k', '>');
    scatter(group.v_prevresp_prct, ones(1,5) * ylm, ...
        sz, 'k', '^');
    axis tight;
  end
end

print(gcf, '-depsc', sprintf('~/Data/serialHDDM/fig3_RTmodulation_correlation.eps'));
disp('printed')
