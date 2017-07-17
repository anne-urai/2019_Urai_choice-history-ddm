function e3_serialBias_SfN_RTmodulation

addpath(genpath('~/code/Tools'));
warning off; close all; clear;
global datasets datasetnames

set(groot, 'defaultaxesfontsize', 6, 'defaultaxestitlefontsizemultiplier', 1, ...
    'defaultaxestitlefontweight', 'bold', ...
    'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters', ...
    'DefaultAxesBox', 'off', ...
    'DefaultAxesTickLength', [0.02 0.05], 'defaultaxestickdir', 'out', 'DefaultAxesTickDirMode', 'manual', ...
    'defaultfigurecolormap', [1 1 1], 'defaultTextInterpreter','tex');
usr = getenv('USER');

% ========================================== %
% MODULATION OF SERIAL CHOICE BIAS
% ========================================== %

plots = {'neutral', 'biased'};

for p = 1:2,
    
    switch p
        case 1
            switch usr
                case 'anne' % local
                    datasets = {'RT_RDK', 'Anke_2afc_neutral', 'NatComm', 'projects/0/neurodec/Data/MEG-PL'};
                case 'aeurai' % lisa/cartesius
                    datasets = {'RT_RDK', 'Anke_neutral', 'NatComm', 'MEG'};
            end
            datasetnames = { {'2AFC, RT', 'n = 22'}, {'2AFC, Braun et al. 2017 bioRxiv', 'n = 22'}, ...
                {'2IFC, Urai et al. 2017 NatComm', 'n = 27'}, {'2IFC, replication', 'n = 61'}};
            
        case 2
            switch usr
                case 'anne' % local
                    datasets = {'Anke_2afc_alternating', 'Anke_2afc_neutral', 'Anke_2afc_repetitive'};
                case 'aeurai' % lisa/cartesius
                    datasets = {'Anke_alternating', 'Anke_neutral', 'Anke_repetitive', 'Anke_serial'};
            end
            datasetnames = {'2AFC alternating', '2AFC neutral', '2AFC repetitive', '2AFC all'};
            datasetnames =  {{'2AFC, Braun et al. 2017 bioRxiv', 'Alternating'}, ...
                {'2AFC, Braun et al. 2017 bioRxiv', 'Neutral'}, ...w
                {'2AFC, Braun et al. 2017 bioRxiv', 'Repetitive'}, ...
                {'2AFC, Braun et al. 2017 bioRxiv', 'All'}};
    end
    
    cnt = 1; close all;
    params = {'dc', 'z'};
    
    for pa = 1:length(params),
        for d = 1:length(datasets),
            
            % get traces for the model with pupil and rt modulation
            traces = readtable(sprintf('~/Data/%s/HDDM/regress_dc_z_prevresp_prevstim_prevrt_prevpupil/all_traces.csv', datasets{d}));
            
            switch params
                case 'dc'
                    dat1 = traces.v_prevresp_prevrt;
                    dat2 = traces.v_prevresp_prevpupil;
                case 'z'
                    dat1 = invlogit(traces.z_prevresp_prevrt);
                    dat2 = invlogit(traces.z_prevresp_prevpupil);
            end
            
            colors = cbrewer('qual', 'Paired', 8);
            % plot the pupil and RT traces
            subplot(4,4,d); hold on;
            h1 = histogram_smooth(dat1, colors(1, :), colors(2, :));
            h2 = histogram_smooth(dat, colors(3, :), colors(4, :));
            
            % show if these are significant (1-sided?)
            pvalRT    = min([mean(dat1 > 0) mean(dat1 < 0)]);
            pvalPupil    = min([mean(dat2 > 0) mean(dat2 < 0)]);
            
            axis tight; axis square;
            l = legend([h1 h2], {sprintf('RT, p = %.3f', pvalRT); sprintf('Pupil, p = %.3f', pvalPupil)}, ...
                'location', 'southeast');
            l.Position(2) = l.Position(2) - 0.15;
            legend boxoff;
            title(datasetnames{d}); xlabel(sprinft('%s ~ prevresp modulation', params{pa}));
            vline(0);
            offsetAxes_y;
        end
    end
    print(gcf, '-depsc', sprintf('~/Data/serialHDDM/figure2_modulationTraces_%s.eps', plots{p}));
    
end
end

function h = histogram_smooth(x, color1, color2)
[f,xi] = ksdensity(x);
area(xi, f, 'edgecolor', 'none', 'facecolor', color1);
% area
h = plot(xi, f, 'color', color2, 'linewidth', 1);
end

function offsetAxes_y()

if ~exist('ax', 'var'), ax = gca;
end
if ~exist('offset', 'var'), offset = 4;
end

% ax.YLim(1) = ax.YLim(1)-(ax.YTick(2)-ax.YTick(1))/offset;
ax.XLim(1) = ax.XLim(1)-(ax.XTick(2)-ax.XTick(1))/offset;

% this will keep the changes constant even when resizing axes
addlistener(ax, 'MarkedClean', @(obj,event)resetVertex(ax));
end

function resetVertex ( ax )
% repeat for Y (set 2nd row)
ax.YRuler.Axle.VertexData(2,1) = min(get(ax, 'Ytick'));
ax.YRuler.Axle.VertexData(2,2) = max(get(ax, 'Ytick'));
% X, Y and Z row of the start and end of the individual axle.
ax.XRuler.Axle.VertexData(1,1) = min(get(ax, 'Xtick'));
ax.XRuler.Axle.VertexData(1,2) = max(get(ax, 'Xtick'));
end
