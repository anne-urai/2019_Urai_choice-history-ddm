function e3_serialBias_SfN_Posteriors_StartingPoint

% Code to fit the history-dependent drift diffusion models described in
% Urai AE, Gee JW de, Donner TH (2018) Choice history biases subsequent evidence accumulation. bioRxiv:251595
%
% MIT License
% Copyright (c) Anne Urai, 2018
% anne.urai@gmail.com

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
	fprintf('%s %s, mean zbias %.4f, range %.4f to %.4f, p-value %.4f \n', datasetnames{d}{1}, datasetnames{d}{2}, ...
	mean(zbias), min(zbias), max(zbias), pval);
    
	close all;
	subplot(4,4,1); hold on;
	h2 = histogram_smooth(zbias, colors(1, :));
    
	% show if these are significant - two sided
	% https://github.com/jwdegee/2017_eLife/blob/master/hddm_regression.py, line 273
    
	axis square; axis tight;
	set(gca, 'xtick', [min(get(gca, 'xlim')) max(get(gca, 'xlim'))], ...
	'xticklabel', []);
	ylim([-0.15 0.05]); set(gca, 'ytick', [-0.15:0.1:0.05]);  
	hline(0);
	% offsetAxes_y;
	offsetAxes;
   
	title(datasetnames{d});
	xlabel('  ');
	txt = sprintf('p = %.4f', pval);
	if pval < 0.0001, txt = 'p < 0.0001'; 
	end
	text(0.1*mean(get(gca, 'xlim')), median(zbias),  ...
	txt, 'horizontalalignment', 'left', 'fontsize', 4);
	set(gca, 'xcolor', 'k', 'ycolor', 'k');
	ylabel('History shift in z', 'color', colors(1, :));

	tightfig;
	print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/zbias_posteriors_d%d.pdf', d));
    
	traces = readtable(sprintf('%s/%s/stimcoding_dc_z_prevresp/group_traces.csv', mypath, datasets{d}));
	zbias = traces.dc_1_ - traces.dc__1_;
	pval = mean(zbias < 0);
	fprintf('%s %s, mean zbias %.4f, range %.4f to %.4f, p-value %.4f \n', datasetnames{d}{1}, datasetnames{d}{2}, ...
	mean(zbias), min(zbias), max(zbias), pval);
    
	close all;
	subplot(4,4,1); hold on;
	h2 = histogram_smooth(zbias, colors(2, :));
    
	% show if these are significant - two sided
	% https://github.com/jwdegee/2017_eLife/blob/master/hddm_regression.py, line 273
    
	axis square; axis tight;
	set(gca, 'xtick', [min(get(gca, 'xlim')) max(get(gca, 'xlim'))], ...
	'xticklabel', []);
	ylim([-0.4 0.8]); set(gca, 'ytick', [-0.4:0.4:0.8]);
	hline(0);
	offsetAxes;
	
	title(datasetnames{d});
	set(gca, 'xcolor', 'k', 'ycolor', 'k');

	ylabel('History shift in v_{bias}', 'color', colors(2, :));
	xlabel('Probability');
	txt = sprintf('p = %.4f', pval);
	if pval < 0.0001, txt = 'p < 0.0001'; 
	end
	text(0.1*mean(get(gca, 'xlim')), median(zbias),  ...
	txt, 'horizontalalignment', 'left', 'fontsize', 4);
	tightfig;
	print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/vbias_posteriors_d%d.pdf', d));
    
end

end

function h = histogram_smooth(x, color2)

[f,xi] = ksdensity(x);
a1 = area(f, xi, 'edgecolor', 'none', 'facecolor', ...
color2, 'facealpha', 0.4, 'showbaseline', 'off', 'basevalue', median(x));

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
