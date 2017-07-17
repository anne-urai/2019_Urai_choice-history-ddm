function e3_serialBias_SfN_Posteriors

addpath(genpath('~/code/Tools'));
warning off; close all; clear;
global datasets datasetnames

set(groot, 'defaultaxesfontsize', 5, 'defaultaxestitlefontsizemultiplier', 1, ...
    'defaultaxestitlefontweight', 'bold', ...
    'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters', ...
    'DefaultAxesBox', 'off', ...
    'DefaultAxesTickLength', [0.02 0.05], 'defaultaxestickdir', 'out', 'DefaultAxesTickDirMode', 'manual', ...
    'defaultfigurecolormap', [1 1 1], 'defaultTextInterpreter','tex');
usr = getenv('USER');

% ========================================== %
% MODULATION OF SERIAL CHOICE BIAS
% use the STIMCODING POSTERIORS!
% ========================================== %

plots = {'neutral', 'biased'};

for p = 2:length(plots),
    
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
    parameters = {'dc', 'z'}; cnt = 1;
    for pa = 1:length(parameters),
        for d = 1:length(datasets),
            
            % get traces for the model with pupil and rt modulation
            traces = readtable(sprintf('~/Data/%s/HDDM/stimcoding_dc_z_prevresp/all_traces.csv', datasets{d}));
            
            switch parameters{pa}
                case 'dc'
                    dat1 = traces.dc_1_;
                    dat2 = traces.dc__1_;
                case 'z' % pull through the inverse logit
                    % see https://groups.google.com/forum/#!topic/hddm-users/PjH7PWZKbqo
                    dat1 = invlogit(traces.z_trans_1_);
                    dat2 = invlogit(traces.z_trans__1_);
            end
            
            % plot the pupil and RT traces
            purples  = cbrewer('seq', 'Oranges', 3);
            greens   = cbrewer('seq', 'Blues', 3);
            
            subplot(4,4,cnt); hold on; cnt = cnt + 1;
            h1 = histogram_smooth(dat1, greens(end, :));
            h2 = histogram_smooth(dat2, purples(end, :));
            
            % show if these are significant - two sided
            % https://github.com/jwdegee/2017_eLife/blob/master/hddm_regression.py, line 273
            pvalD    = min([mean(dat2 > dat1) mean(dat2 < dat1)]);
            
            axis tight; axis square;
            
            txt = sprintf('p = %.3f', pvalD);
            if pvalD < 0.001,
                txt = sprintf('p < 0.001');
            end
            text(min(get(gca, 'xlim')) + 0.75*(range(get(gca, 'xlim'))), ...
                min(get(gca, 'ylim')) + 0.8*(range(get(gca, 'ylim'))), ...
                txt, 'fontsize', 5);
            
            if pa == 1,
                title(datasetnames{d});
            end
            xlabel(parameters{pa});
            offsetAxes_y;
        end
        if cnt == 4, cnt = 5 ; end
        % cnt = cnt + 4;
    end
    print(gcf, '-depsc', sprintf('~/Data/serialHDDM/figure2_posteriors_%s_stimcoding.eps', plots{p}));
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/figure2_posteriors_%s_stimcoding.pdf', plots{p}));
    
end
end

function h = histogram_smooth(x, color2)

[f,xi] = ksdensity(x);
a1 = area(xi, f, 'edgecolor', 'none', 'facecolor', color2, 'facealpha', 0.4);

% % Make area transparent
% drawnow; % pause(0.05);  % This needs to be done for transparency to work
% a1.Face.ColorType = 'truecoloralpha';
% a1.Face.ColorData(4) = 255 * 0.3; % Your alpha value is the 0.3

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
