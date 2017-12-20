function e3_serialBias_SfN_Posteriors_StartingPoint

	addpath(genpath('~/code/Tools'));
	warning off; close all;
	global mypath datasets datasetnames colors

	% ========================================== %
	% MODULATION OF SERIAL CHOICE BIAS
	% use the STIMCODING POSTERIORS!
	% ========================================== %

	for d = 1:length(datasets),
		traces = readtable(sprintf('%s/%s/stimcoding_dc_z_prevresp/group_traces.csv', mypath, datasets{d}));
		zbias = invlogit(traces.z_trans_1_) - invlogit(traces.z_trans__1_);
		pval = mean(zbias > 0);
		fprintf('%s %s, mean zbias %.4f, range %.4f to %.4f, p-value %.4f \n', datasetnames{d}{1}{1}, datasetnames{d}{1}{2}, ...
		 mean(zbias), min(zbias), max(zbias), pval);

    close all;
    subplot(4,4,1); hold on;
    h2 = histogram_smooth(zbias, colors(1, :));
    
    % show if these are significant - two sided
    % https://github.com/jwdegee/2017_eLife/blob/master/hddm_regression.py, line 273
    
    axis square; axis tight;
    set(gca, 'ytick', [min(get(gca, 'ytick')) max(get(gca, 'ytick'))], ...
        'yticklabel', []);
		if max(get(gca, 'xlim')) < 0, xlim([min(get(gca, 'xlim')) 0]);
		disp(get(gca, 'xlim')); end
	set(gca, 'xtick', unique(sort([min(get(gca, 'xtick')) 0 max(get(gca, 'xtick'))])));
    offsetAxes_y;
    title(datasetnames{d}{1});
	xlabel('History shift in z');
	ylabel('Probability');
	txt = sprintf('p = %.4f', pval);
	if pval < 0.0001, txt = 'p < 0.0001'; end
	text(0.4*median(zbias), 0.1*mean(get(gca, 'ylim')), ...
	txt, 'horizontalalignment', 'center', 'fontsize', 4);
    tightfig;
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/zbias_posteriors_d%d.pdf', d));
	
end

end

function h = histogram_smooth(x, color2)

[f,xi] = ksdensity(x);
a1 = area(xi, f, 'edgecolor', 'none', 'facecolor', ...
    color2, 'facealpha', 0.4, 'showbaseline', 'off');

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
