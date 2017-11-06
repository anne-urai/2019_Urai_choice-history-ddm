function e3_serialBias_SfN_RTmodulation

addpath(genpath('~/code/Tools'));
warning off; close all; clear;
global datasets datasetnames

% ========================================== %
% MODULATION OF SERIAL CHOICE BIAS
% ========================================== %

close all;
for d = 1:length(datasets),
    
    % get traces for the model with pupil and rt modulation
    traces = readtable(sprintf('~/Data/HDDM/%s/regress_dc_z_prevresp_prevrt/all_traces.csv', datasets{d}));
    
    dat1 = traces.z_prevresp_prevrt;
    dat2 = traces.v_prevresp_prevrt;
 
    % colors
    colors = [141 165 8;  8 141 165; 150 150 150] ./ 256;
    
    % plot the pupil and RT traces
    subplot(4,4,d); hold on;
    h1 = histogram_smooth(dat1, [], colors(1, :));
    h2 = histogram_smooth(dat2, [], colors(2, :));
    
    % show if these are significant (1-sided?)
    pvalZ       = min([mean(dat1 > 0) mean(dat1 < 0)]);
    pvalV       = min([mean(dat2 > 0) mean(dat2 < 0)]);
    
    axis tight; axis square;
    l = legend([h1 h2], {sprintf('zbias, p = %.3f', pvalZ); sprintf('vbias, p = %.3f', pvalV)}, ...
        'location', 'southeast');
    l.Position(2) = l.Position(2) - 0.15;
    legend boxoff;
    title(datasetnames{d}); 
    %xlabel(sprinft('%s ~ prevresp modulation', params{pa}));
    %vline(0);
    %offsetAxes_y;
    print(gcf, '-depsc', sprintf('~/Data/serialHDDM/figure2_modulationTraces_D%s.eps', d));
    
end
end

function h = histogram_smooth(x, color1, color2)
[f,xi] = ksdensity(x);
if ~isempty(color1),
    area(xi, f, 'edgecolor', 'none', 'facecolor', color1);
end
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
