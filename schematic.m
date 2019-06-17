function schematic(seed)

% Code to fit the history-dependent drift diffusion models as described in
% Urai AE, de Gee JW, Tsetsos K, Donner TH (2019) Choice history biases subsequent evidence accumulation. eLife, in press.
%
% MIT License
% Copyright (c) Anne Urai, 2019
% anne.urai@gmail.com

close all; clc;

global colors
% colors = [8 141 165; 141 165 8;  150 150 150] ./ 256;
pairedcols = cbrewer('qual', 'Paired', 8);

fz = 5; timefz = 5;
set(groot, 'defaultaxesfontsize', fz, 'defaultaxestitlefontsizemultiplier', 1);

%%  random seed will determine how the drifting timecourse looks
if ~exist('seed', 'var'), seed = 100:120; end
disp(seed);

% use code by Peter Murphy to compute RT distributions
addpath('analyticalDDM/DDM');
addpath(genpath('~/code/Tools'));

tmax = 0.7;
%% ===================================== %%
% make 3 sets of RT distributions
% without any bias
% ===================================== %%

pm = [0.02 0 0.1 0.001 0.05]; % standard parameters
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
cfg.timestep = 0.005; % 100 ms
cfg.time     = cfg.timestep:cfg.timestep:tmax;
cfg.a        = 1; % bound with, z = 0
cfg.cdW      = 0.1; % variance of normally distributed noise
cfg.v        = 2.5*cfg.timestep; % drift rate per timestep
cfg.z        = 0; % drift rate per timestep
defcfg = cfg;

%% ===================================== %%
% make an overview of the two biasing mechanisms in the DDM
% ===================================== %%

close all; subplot(431); hold on;
for s = 1:length(seed)
    cfg.seed     = s;
    y1 = ddm(cfg);
    % plot the drift on top
    plot(cfg.time, y1, 'color', [0.5 0.5 0.5], 'linewidth', 0.2);
end

cfg.v = cfg.v+1*cfg.timestep;
for s = 1:length(seed)
    cfg.seed     = s+100;
    y1 = ddm(cfg);
    % plot the drift on top
    plot(cfg.time, y1, 'color', pairedcols(1, :), 'linewidth', 0.2);
end

% show the unbiased average drift towards two stimuli
cfg.v = 2.5*cfg.timestep;
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
plot(cfg.time, y, 'color', pairedcols(2, :));
cfg.v = -cfg.v + 2*cfg.dc; % flip around drift rate
y = ddm(cfg);
plot(cfg.time, y,'k', 'color', pairedcols(2, :));

% add distributions at the top!
scaling = 300;
plot(ts, scaling*gC_nobias + cfg.a, 'k');
meanmedianmode(ts, scaling*gC_nobias + cfg.a, 'k');
plot(ts, scaling*gC_dc + cfg.a, 'color', pairedcols(2, :));
meanmedianmode(ts, scaling*gC_dc + cfg.a, pairedcols(2, :));

plot(ts, -scaling*gE_nobias - cfg.a, 'k');
meanmedianmode(ts, -scaling*gE_nobias - cfg.a, 'k');
plot(ts, -scaling*gE_dc - cfg.a, 'color', pairedcols(2, :));
meanmedianmode(ts, -scaling*gE_dc - cfg.a, pairedcols(2, :));

%%  layout
%ylim([-cfg.a cfg.a]);
axis tight;
set(gca, 'ytick', [-cfg.a cfg.z cfg.a], 'yticklabel', {'0', 'z', 'a'});
text(0.83*max(cfg.time), -0.2, 'Time', 'fontsize', timefz);
%title({'History shift' 'in drift bias'}, 'color', colors(2, :));
%title({'Drift bias (v_{bias})'}, 'color', colors(2, :));
arrow([cfg.time(1) cfg.z ], [cfg.time(end) cfg.z], 'linestyle', ':', 'linewidth', 0.5, 'length', 4, 'TipAngle', 45);

% add two axes manually
plot([cfg.time(1) cfg.time(end)], [cfg.a cfg.a], 'k', 'linewidth', 0.5);
plot([cfg.time(1) cfg.time(end)], [-cfg.a -cfg.a], 'k', 'linewidth', 0.5);
box off;
ax = gca;
addlistener(ax, 'MarkedClean', @(obj,event)resetVertex(ax));
set(ax, 'xcolor', 'w', 'xtick',[]);
xlim([min(cfg.time) max(cfg.time)]);
tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/DDMschematic_driftbias.pdf'));

%% ===================================== %%
% now change in starting point
% ===================================== %%

cfg = defcfg;
close all; subplot(431); hold on;

for s = 1:length(seed)
    cfg.seed     = s;
    y1 = ddm(cfg);
    % plot the drift on top
    plot(cfg.time, y1, 'color', [0.5 0.5 0.5], 'linewidth', 0.2);
end

cfg.z = 0.4;
for s = 1:length(seed)
    cfg.seed     = s+100;
    y1 = ddm(cfg);
    % plot the drift on top
    plot(cfg.time, y1, 'color', pairedcols(3, :), 'linewidth', 0.2);
end

% show the unbiased average drift towards two stimuli
cfg.cdW = 0;
cfg.z = 0;
y = ddm(cfg);
plot(cfg.time, y,'k');
cfg.v = -cfg.v; % flip around drift rate
y = ddm(cfg);
plot(cfg.time, y,'k');

% now with drift criterion bias
cfg.z = 0.4;
y = ddm(cfg);
plot(cfg.time, y, 'color', pairedcols(4, :));
cfg.v = -cfg.v;
y = ddm(cfg);
plot(cfg.time, y,'k', 'color', pairedcols(4, :));

% drifting particle
%plot(cfg.time, y1, 'color', [0.5 0.5 0.5]);

% add distributions at the top!
%plot(ts, scaling*gC_dc + cfg.a, 'color', colors(2, :), 'linestyle', '-.');
%plot(ts, -scaling*gE_dc - cfg.a, 'color', colors(2, :), 'linestyle', '-.');

plot(ts, scaling*gC_nobias + cfg.a, 'k');
meanmedianmode(ts, scaling*gC_nobias + cfg.a, 'k');

plot(ts, scaling*gC_z + cfg.a, 'color', pairedcols(4, :));
meanmedianmode(ts, scaling*gC_z + cfg.a, pairedcols(4, :));

plot(ts, -scaling*gE_nobias - cfg.a, 'k');
meanmedianmode(ts, -scaling*gE_nobias - cfg.a, 'k');

plot(ts, -scaling*gE_z - cfg.a, 'color', pairedcols(4, :));
meanmedianmode(ts, -scaling*gE_z - cfg.a, pairedcols(4, :));


% layout
%ylim([-cfg.a cfg.a]);
axis tight;
set(gca, 'ytick', [-cfg.a 0 cfg.a], 'yticklabel', {'0', 'z', 'a'});
text(0.83*max(cfg.time), -0.2, 'Time', 'fontsize', timefz);
%title({'History shift' 'in starting point'}, 'color', colors(1, :));
plot([cfg.time(1) cfg.time(end)], [cfg.a cfg.a], 'k', 'linewidth', 0.5);
plot([cfg.time(1) cfg.time(end)], [-cfg.a -cfg.a], 'k', 'linewidth', 0.5);
%title({'Starting point (z)'}, 'color', colors(1, :));
arrow([cfg.time(1) 0], [cfg.time(end) 0], 'linestyle', ':', 'linewidth', 0.5, 'length', 4, 'TipAngle', 45);

box off;
ax = gca;
addlistener(ax, 'MarkedClean', @(obj,event)resetVertex(ax));
set(ax, 'xcolor', 'w', 'xtick',[]);
xlim([min(cfg.time) max(cfg.time)]);

tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/DDMschematic_startingpoint.pdf'));

%% now add the equations!
%
% subplot(335);
% xoffset = 0.02;
% text(xoffset, 1.15, 'dy = (s\cdotv+\bf{v_{bias}}\rm{)dt + cdW}', 'fontsize', fz);
% text(xoffset, 1, 'y(0) = z = a/2', 'fontsize', fz);
% axis off;
%
% subplot(334);
% text(xoffset, 1.15, 'dy = s\cdotv\cdotdt + cdW', 'fontsize', fz);
% text(xoffset, 1, 'y(0) = z = a/2 + \bf{z_{bias}}', 'fontsize', fz);
% axis off;

% save
% offsetAxes;


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
firstPassage = min([find(y > cfg.a, 1, 'first') find(y < -cfg.a, 1, 'first')]);
y(firstPassage+1:end) = NaN;
%y(y > cfg.a) = NaN;
%y(y < -cfg.a) = NaN;

end


function resetVertex ( ax )
% extract the x axis vertext data

% repeat for Y (set 2nd row)
ax.YRuler.Axle.VertexData(2,1) = min(get(ax, 'Ytick'));
ax.YRuler.Axle.VertexData(2,2) = max(get(ax, 'Ytick'));

end


function meanmedianmode(x,y, col)

% from cumulative distribution function compute mean, median and mode
% 
% plot(x(dsearchn(y', mean(y))), y(dsearchn(y', mean(y))), 'o', 'markerfacecolor', col, 'markeredgecolor', 'w', 'markersize', 4);
% plot(x(dsearchn(y', median(y))), y(dsearchn(y', median(y))), 's', 'markerfacecolor', col, 'markeredgecolor', 'w', 'markersize', 4);

% mode
if mean(y) > 0,
    plot(x(dsearchn(y', max(y))), y(dsearchn(y', max(y))), 'o', 'markerfacecolor', col, 'markeredgecolor', col, 'markersize', 1);

else
    plot(x(dsearchn(y', min(y))), y(dsearchn(y', min(y))), 'o', 'markerfacecolor', col, 'markeredgecolor', col, 'markersize', 1);
end

end

