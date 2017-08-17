function e3_serialBias_SfN_Posteriors

addpath(genpath('~/code/Tools'));
warning off; close all; clear;
global mypath

% neutral vs biased plots
datasets = {'Anke_2afc_alternating', ...
    'Anke_2afc_neutral', 'Anke_2afc_repetitive'};
datasetnames = {{'2AFC, Braun et al. 2017', 'Alternating'}, ...
    {'2AFC, Braun et al. 2017', 'Neutral'}, ...
    {'2AFC, Braun et al. 2017', 'Repetitive'}};

% ========================================== %
% MODULATION OF SERIAL CHOICE BIAS
% use the STIMCODING POSTERIORS!
% ========================================== %

parameters = {'dc', 'z'};
for pa = 1:length(parameters),
    for d = 1:length(datasets),
        disp(datasets{d});
        % get traces for the model with pupil and rt modulation
        traces = readtable(sprintf('%s/%s/stimcoding_dc_z_prevresp/all_traces.csv', mypath, datasets{d}));
        
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
        colors = linspecer(9, 'sequential');
        
        close all;
        subplot(4,4,1); hold on;
        h1 = histogram_smooth(dat1, colors(3, :));
        h2 = histogram_smooth(dat2, colors(2, :));
        
        % show if these are significant - two sided
        % https://github.com/jwdegee/2017_eLife/blob/master/hddm_regression.py, line 273
        pvalD    = min([mean(dat2 > dat1) mean(dat2 < dat1)]);
        
        axis tight; axis square;
        ylims = get(gca, 'ylim');
        ylim([ylims(1) ylims(2)*1.2]);
        
        offsetAxes_y;
        
        switch parameters{pa}
            case 'dc'
                xlim([-1 1]);
                xlabel('v_{bias}');
                
            case 'z'
                xlim([0.4 0.6]);
                xlabel('z_{bias}');
                
        end
        
        txt = sprintf('p = %.3f', pvalD);
        if pvalD < 0.001,
            txt = sprintf('p < 0.001');
        end
        text(min(get(gca, 'xlim')) + 0.7*(range(get(gca, 'xlim'))), ...
            min(get(gca, 'ylim')) + 0.8*(range(get(gca, 'ylim'))), ...
            txt, 'fontsize', 5);
        
        if pa == 1,
            title(datasetnames{d});
        end
        % offsetAxes_y;
        if d == 1, ylabel('Posterior probability'); end
        ticks = get(gca, 'ytick'); set(gca, 'ytick', ticks([1 end]));
        tightfig;
        
        print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/figure2_posteriors_%s_d%d.pdf', parameters{pa}, d));
    end
end
close all;

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
set(gca, 'color', 'none');

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
