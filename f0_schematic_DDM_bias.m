function f0_schematic_DDM_bias(seed)
close all; clc;

colors = [8 141 165; 141 165 8;  150 150 150] ./ 256;

fz = 9; timefz = 8;
set(groot, 'defaultaxesfontsize', fz, 'defaultaxestitlefontsizemultiplier', 1);

% random seed will determine how the drifting timecourse looks
if ~exist('seed', 'var'), seed = 100; end
disp(seed);

% use code by Peter Murphy to compute RT distributions
addpath('/Users/anne/Drive/Dropbox/code/analyticalDDM/DDM');

tmax = 0.7;
%% make 3 sets of distributions
% without any bias
pm = [0.02 0 0.1 0.001 0.05];
[gC_nobias,gE_nobias,ts] = fpt_regular_DDM(pm, tmax);
plot(ts, gC_nobias, ts, gE_nobias);

% biased drift towards option 1
pm(1) = pm(1) * 10;
[gC_dc,gE_dc,ts] = fpt_regular_DDM(pm, tmax);
plot(ts, gC_dc, ts, gE_dc);

% biased starting point towards option 1
pm = [0.02 0 0.1 0.001 0.05];
pm(5) = pm(5) * 1.4;
[gC_z,gE_z,ts] = fpt_regular_DDM(pm, tmax);
plot(ts, gC_z, ts, gE_z);

%% set parameters
cfg.timestep = 0.01; % 100 ms
cfg.time     = cfg.timestep:cfg.timestep:tmax;
cfg.a        = 1; % bound with, z = 0
cfg.cdW      = 0.1; % variance of normally distributed noise
cfg.v        = 2.5*cfg.timestep; % drift rate per timestep
cfg.z        = 0; % drift rate per timestep
cfg.seed     = seed;
defcfg = cfg;

%% make an overview of the two biasing mechanisms in the DDM

subplot(331); hold on;
arrow([cfg.time(1) cfg.z ], [cfg.time(end) cfg.z], 'linewidth', 0.5, 'length', 4, 'TipAngle', 45);
y1 = ddm(cfg);

% show the unbiased average drift towards two stimuli
cfg.cdW = 0;
y = ddm(cfg);
plot(cfg.time, y,'k');
cfg.v = -cfg.v; % flip around drift rate
y = ddm(cfg);
plot(cfg.time, y,'k');

% now with drift criterion bias
cfg.dc = 1*cfg.timestep;
cfg.v = cfg.v+cfg.dc;
y = ddm(cfg);
plot(cfg.time, y, 'color', colors(1, :));
cfg.v = -cfg.v + 2*cfg.dc; % flip around drift rate
y = ddm(cfg);
plot(cfg.time, y,'k', 'color', colors(1, :));

% plot the drift on top 
plot(cfg.time, y1, 'color', [0.5 0.5 0.5]);

% add distributions at the top!
scaling = 300;
plot(ts, -scaling*gE_z - cfg.a, 'color', colors(2, :), 'linestyle', '-.');
plot(ts, scaling*gC_z + cfg.a, 'color', colors(2, :), 'linestyle', '-.');

plot(ts, scaling*gC_nobias + cfg.a, 'k');
plot(ts, scaling*gC_dc + cfg.a, 'color', colors(1, :));
plot(ts, -scaling*gE_nobias - cfg.a, 'k');
plot(ts, -scaling*gE_dc - cfg.a, 'color', colors(1, :));

%%  layout
%ylim([-cfg.a cfg.a]); 
axis tight;
set(gca, 'ytick', [-cfg.a cfg.z cfg.a], 'yticklabel', {'0', 'z', 'a'});
text(0.83*max(cfg.time), -0.2, 'Time', 'fontsize', timefz);
title('Biased drift');
% add two axes manually
plot([cfg.time(1) cfg.time(end)], [cfg.a cfg.a], 'k', 'linewidth', 0.5);
plot([cfg.time(1) cfg.time(end)], [-cfg.a -cfg.a], 'k', 'linewidth', 0.5);
box off; 
ax = gca;
addlistener(ax, 'MarkedClean', @(obj,event)resetVertex(ax));
set(ax, 'xcolor', 'w');
xlim([min(cfg.time) max(cfg.time)]);

%% now change in starting point
cfg = defcfg;
subplot(332); hold on;
arrow([cfg.time(1) cfg.z ], [cfg.time(end) cfg.z], 'linewidth', 0.5, 'length', 4, 'TipAngle', 45);

% show the unbiased average drift towards two stimuli
cfg.cdW = 0;
y = ddm(cfg);
plot(cfg.time, y,'k');
cfg.v = -cfg.v; % flip around drift rate
y = ddm(cfg);
plot(cfg.time, y,'k');

% now with drift criterion bias
cfg.z = 0.4;
y = ddm(cfg);
plot(cfg.time, y, 'color', colors(2, :));
cfg.v = -cfg.v;
y = ddm(cfg);
plot(cfg.time, y,'k', 'color', colors(2, :));

% drifting particle
plot(cfg.time, y1, 'color', [0.5 0.5 0.5]);

% add distributions at the top!
plot(ts, scaling*gC_dc + cfg.a, 'color', colors(1, :), 'linestyle', '-.');
plot(ts, -scaling*gE_dc - cfg.a, 'color', colors(1, :), 'linestyle', '-.');

plot(ts, scaling*gC_nobias + cfg.a, 'k');
plot(ts, scaling*gC_z + cfg.a, 'color', colors(2, :));
plot(ts, -scaling*gE_nobias - cfg.a, 'k');
plot(ts, -scaling*gE_z - cfg.a, 'color', colors(2, :));

% layout
%ylim([-cfg.a cfg.a]);
axis tight;
set(gca, 'ytick', [-cfg.a 0 cfg.a], 'yticklabel', {'0', 'z', 'a'});
text(0.83*max(cfg.time), -0.2, 'Time', 'fontsize', timefz);
title('Biased starting point');
plot([cfg.time(1) cfg.time(end)], [cfg.a cfg.a], 'k', 'linewidth', 0.5);
plot([cfg.time(1) cfg.time(end)], [-cfg.a -cfg.a], 'k', 'linewidth', 0.5);

box off;
ax = gca;
addlistener(ax, 'MarkedClean', @(obj,event)resetVertex(ax));
set(ax, 'xcolor', 'w');
xlim([min(cfg.time) max(cfg.time)]);

%% now add the equations!

subplot(334);
xoffset = 0.02;
text(xoffset, 1.15, 'dy = (s\cdotv+\bf{v_{bias}}\rm{)dt + cdW}', 'fontsize', fz);
text(xoffset, 1, 'y(0) = z = a/2', 'fontsize', fz);
axis off;

subplot(335);
text(xoffset, 1.15, 'dy = s\cdotv\cdotdt + cdW', 'fontsize', fz);
text(xoffset, 1, 'y(0) = z = a/2 + \bf{z_{bias}}', 'fontsize', fz);
axis off;

% save
% offsetAxes;
tightfig;
print(gcf, '-depsc', sprintf('~/Data/serialHDDM/DDMschematic.eps'));
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/DDMschematic.pdf'));

end


function y = ddm(cfg)

% generate a sample of evidence for each moment in time
evidence = cfg.v * ones(size(cfg.time));

% add noise
rng(cfg.seed); % make sure this is the same each time
noise = randn(size(evidence)) .* cfg.cdW;
noise(1) = cfg.z; % start with no noise
noisyevidence = evidence + noise;

integratedevidence = cumsum(noisyevidence);
y = integratedevidence;
firstPassage = find(y > cfg.a);
y(firstPassage:end) = NaN;
y(y > cfg.a) = NaN;
y(y < -cfg.a) = NaN;

end


function resetVertex ( ax )
% extract the x axis vertext data

% repeat for Y (set 2nd row)
ax.YRuler.Axle.VertexData(2,1) = min(get(ax, 'Ytick'));
ax.YRuler.Axle.VertexData(2,2) = max(get(ax, 'Ytick'));

end

%% additional plots
