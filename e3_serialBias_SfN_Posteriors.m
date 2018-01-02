function e3_serialBias_SfN_Posteriors

addpath(genpath('~/code/Tools'));
warning off; close all;
global mypath datasets datasetnames

% neutral vs biased plots
% datasets = {'Anke_2afc_alternating', ...
%     'Anke_2afc_neutral', 'Anke_2afc_repetitive'};
% datasetnames = {{'2AFC, Braun et al. 2017', 'Alternating'}, ...
%     {'2AFC, Braun et al. 2017', 'Neutral'}, ...
%     {'2AFC, Braun et al. 2017', 'Repetitive'}};

datasets =     {'Anke_2afc_sequential', 'Anke_MEG', 'Bharath_fMRI'};

datasetnames = {   {'Anke JoN'}, ...
   {'Anke MEG'}, {'Bharath fMRI'}};

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
                dat1 = traces.dc_1_0_20_0_ - traces.dc__1_0_20_0_;
                dat2 = traces.dc_1_0_50_0_ - traces.dc__1_0_50_0_;
                dat3 = traces.dc_1_0_80_0_ - traces.dc__1_0_80_0_;
                
            case 'z' % pull through the inverse logit
                % see https://groups.google.com/forum/#!topic/hddm-users/PjH7PWZKbqo
                dat1 = invlogit(traces.z_trans_1_0_20_0_) - invlogit(traces.z_trans__1_0_20_0_);
                dat2 = invlogit(traces.z_trans_1_0_50_0_) - invlogit(traces.z_trans__1_0_50_0_);
                dat3 = invlogit(traces.z_trans_1_0_80_0_) - invlogit(traces.z_trans__1_0_80_0_);
        end
        
        % color in different grouos
        colors = cbrewer('qual', 'Paired', 10);
        transitioncolors = [[0.5 0.5 0.5]; colors([7 9], :)];
        meancolors = [0 0 0; colors([8 10], :)];
        
        close all;
        subplot(4,4,1); hold on;
        h2 = histogram_smooth(dat1, transitioncolors(2, :));
        h1 = histogram_smooth(dat2, transitioncolors(1, :));
        h2 = histogram_smooth(dat3, transitioncolors(3, :));
        
        % show if these are significant - two sided
        % https://github.com/jwdegee/2017_eLife/blob/master/hddm_regression.py, line 273
        
        axis tight; axis square;
        xlims = get(gca, 'xlim');
        xlim([xlims(1) xlims(2)*1.2]);
        set(gca, 'xtick', [0 max(get(gca, 'xtick'))]);
        offsetAxes_y;
        
        switch parameters{pa}
            case 'dc'
               % ylim([-1 1]);
                if d == 1, ylabel('History shift in v'); end
                set(gca, 'ytick', [-1 0 1], 'ylim', [-1.3 1.3]);
            case 'z'
                % ylim([-0.1 0.1]);
                if d == 1, ylabel('History shift in z'); end
                set(gca, 'ytick', [-.1 0 .1], 'ylim', [-0.12 0.12]);
        end
        
        pvalD = posteriorpval(dat1, dat2);
        txt = sprintf('p = %.3f', pvalD);
        if pvalD < 0.001,
            txt = sprintf('p < 0.001');
        end
        text(min(get(gca, 'xlim')) + 0.05*(range(get(gca, 'xlim'))), ...
            nanmean(dat1), txt, 'fontsize', 5, 'color', meancolors(2, :));
        
        pvalD = posteriorpval(dat3, dat2);
        txt = sprintf('p = %.3f', pvalD);
        if pvalD < 0.001,
            txt = sprintf('p < 0.001');
        end
        text(min(get(gca, 'xlim')) + 0.05*(range(get(gca, 'xlim'))), ...
           nanmean(dat3), txt, 'fontsize', 5, 'color', meancolors(3, :));
       
        if pa == 1,
            title(datasetnames{d});
        end
        
        if pa == 2, xlabel('Posterior probability'); end
        ticks = get(gca, 'xtick'); set(gca, 'xtick', [0 ticks(end)]);
        set(gca, 'xlim', [0 max(get(gca, 'xlim'))]);
        offsetAxes_y;

        tightfig;
        print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/figure2_posteriors_%s_d%d.pdf', parameters{pa}, d));
  
    end
end
close all;

end

function h = histogram_smooth(x, color2)

[f,xi] = ksdensity(x);
% a1 = area(f, xi, 'edgecolor', 'none', 'facecolor', ...
%     color2, 'facealpha', 0.4, 'showbaseline', 'off', 'basevalue', median(xi));
patch(f,xi, color2, 'edgecolor', 'none', 'facealpha', 0.4);

% % Make area transparent
% drawnow; % pause(0.05);  % This needs to be done for transparency to work
% a1.Face.ColorType = 'truecoloralpha';
% a1.Face.ColorData(4) = 255 * 0.3; % Your alpha value is the 0.3

% area
h = plot(f, xi, 'color', color2, 'linewidth', 1);
set(gca, 'color', 'none');

end

function offsetAxes_y()

if ~exist('ax', 'var'), ax = gca;
end
if ~exist('offset', 'var'), offset = 4;
end

% ax.YLim(1) = ax.YLim(1)-(ax.YTick(2)-ax.YTick(1))/offset;
ax.YLim(2) = ax.YLim(2)+(ax.YTick(2)-ax.YTick(1))/offset;

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
