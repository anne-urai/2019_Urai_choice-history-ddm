function schematic(seed)

% Code to fit the history-dependent drift diffusion models described in
% Urai AE, Gee JW de, Donner TH (2018) Choice history biases subsequent evidence accumulation. bioRxiv:251595
%
% MIT License
% Copyright (c) Anne Urai, 2018
% anne.urai@gmail.com

close all; clc;

global colors
% colors = [8 141 165; 141 165 8;  150 150 150] ./ 256;

fz = 8; timefz = 6;
set(groot, 'defaultaxesfontsize', fz, 'defaultaxestitlefontsizemultiplier', 1);

%% ===================================== %%
% CONDITIONAL BIAS FUNCTIONS
% ===================================== %%

% close all;
% subplot(441); hold on;
% CBF(0, 0.5, [0.3 0.3 0.3]); % no history bias
% CBF(0.05, 0.5, colors(2, :)); % starting point bias
% CBF(0, 0.625, colors(1, :)); % drift bias
% 
% % layout
% set(gca, 'xticklabel', [], 'fontsize', timefz);
% ylabel('Choice bias');
% xlabel({'Response time'});
% set(gca, 'yticklabel', []);
% axis tight; offsetAxes; tightfig;
% print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/CBF_schematic.pdf'));
% close all;

%%  random seed will determine how the drifting timecourse looks
if ~exist('seed', 'var'), seed = 100; end
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
cfg.timestep = 0.01; % 100 ms
cfg.time     = cfg.timestep:cfg.timestep:tmax;
cfg.a        = 1; % bound with, z = 0
cfg.cdW      = 0.1; % variance of normally distributed noise
cfg.v        = 2.5*cfg.timestep; % drift rate per timestep
cfg.z        = 0; % drift rate per timestep
cfg.seed     = seed;
defcfg = cfg;

%% ===================================== %%
% make an overview of the two biasing mechanisms in the DDM
% ===================================== %%

close all; subplot(441); hold on;
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
plot(cfg.time, y, 'color', colors(2, :));
cfg.v = -cfg.v + 2*cfg.dc; % flip around drift rate
y = ddm(cfg);
plot(cfg.time, y,'k', 'color', colors(2, :));

% plot the drift on top
plot(cfg.time, y1, 'color', [0.5 0.5 0.5]);

% add distributions at the top!
scaling = 300;
%plot(ts, -scaling*gE_z - cfg.a, 'color', colors(1, :), 'linestyle', '-.');
%plot(ts, scaling*gC_z + cfg.a, 'color', colors(1, :), 'linestyle', '-.');

plot(ts, scaling*gC_nobias + cfg.a, 'k');
plot(ts, scaling*gC_dc + cfg.a, 'color', colors(2, :));
plot(ts, -scaling*gE_nobias - cfg.a, 'k');
plot(ts, -scaling*gE_dc - cfg.a, 'color', colors(2, :));

%%  layout
%ylim([-cfg.a cfg.a]);
axis tight;
set(gca, 'ytick', [-cfg.a cfg.z cfg.a], 'yticklabel', {'0', 'z', 'a'});
text(0.83*max(cfg.time), -0.2, 'Time', 'fontsize', timefz);
%title({'History shift' 'in drift bias'}, 'color', colors(2, :));
%title({'Drift bias (v_{bias})'}, 'color', colors(2, :));

% add two axes manually
plot([cfg.time(1) cfg.time(end)], [cfg.a cfg.a], 'k--', 'linewidth', 0.5);
plot([cfg.time(1) cfg.time(end)], [-cfg.a -cfg.a], 'k--', 'linewidth', 0.5);
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
close all; subplot(441); hold on;
arrow([cfg.time(1) cfg.z], [cfg.time(end) cfg.z], 'linewidth', 0.5, 'length', 4, 'TipAngle', 45);

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
plot(cfg.time, y, 'color', colors(1, :));
cfg.v = -cfg.v;
y = ddm(cfg);
plot(cfg.time, y,'k', 'color', colors(1, :));

% drifting particle
plot(cfg.time, y1, 'color', [0.5 0.5 0.5]);

% add distributions at the top!
%plot(ts, scaling*gC_dc + cfg.a, 'color', colors(2, :), 'linestyle', '-.');
%plot(ts, -scaling*gE_dc - cfg.a, 'color', colors(2, :), 'linestyle', '-.');

plot(ts, scaling*gC_nobias + cfg.a, 'k');
plot(ts, scaling*gC_z + cfg.a, 'color', colors(1, :));
plot(ts, -scaling*gE_nobias - cfg.a, 'k');
plot(ts, -scaling*gE_z - cfg.a, 'color', colors(1, :));

% layout
%ylim([-cfg.a cfg.a]);
axis tight;
set(gca, 'ytick', [-cfg.a 0 cfg.a], 'yticklabel', {'0', 'z', 'a'});
text(0.83*max(cfg.time), -0.2, 'Time', 'fontsize', timefz);
%title({'History shift' 'in starting point'}, 'color', colors(1, :));
plot([cfg.time(1) cfg.time(end)], [cfg.a cfg.a], 'k--', 'linewidth', 0.5);
plot([cfg.time(1) cfg.time(end)], [-cfg.a -cfg.a], 'k--', 'linewidth', 0.5);
%title({'Starting point (z)'}, 'color', colors(1, :));

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


%% COMPUTE CONDITIONAL BIAS FUNCTIONS %%
function pbias = CBF(dc, z, color)

time = 0.01:0.01:2;

% from https://github.com/anne-urai/serialDDM/blob/master/simulations/1_ddm_rts.py
v = 0.1;
a = 0.1;
% ndt = 0.1;
% dc_slope = 0;
% ou = 0;;
dt = 0.01;
scaling = 0.1;

% for X trials, generate RTs and responses
% generate a sample of evidence for each moment in time
nsim        = 100000;
dv          = nan(nsim, numel(time));
dv(:, 1)    = z*a; % set starting point
for t = 2:length(time),
    dv(:, t) = dv(:, t-1) + (v+dc) .* dt + ...
        normrnd(0, 1, [nsim 1]) .* sqrt(dt) .* scaling;
end

% apply bounds
RT  = nan(1, nsim);
choice = nan(1, nsim);

for n = 1:nsim,
    for j = 1:length(time),
        if dv(n,j) >= a,
            RT(n) = time(j);
            choice(n) = 1;
            break
        elseif dv(n,j) <= 0,
            RT(n) = time(j);
            choice(n) = 0;
            break
        end
    end
end

% for i = 1:nsim,
%     
%     evidence            = (v+dc) .* dt .* ones(size(time));
%     noise               = normrnd(0, 1, size(evidence)) .* sqrt(dt) .* scaling;
%     noisyevidence       = evidence + noise;
%     noisyevidence(1)    = noisyevidence(1) + z*a;
%     integratedevidence  = cumsum(noisyevidence);
%     %integratedevidence  = cumsum(noisyevidence) + z*a;
%     
%     % apply bound
%     [firstPassage, boundPassed] = min([find(integratedevidence >= a, 1) find(integratedevidence <= 0, 1)]);
%     if ~isempty(firstPassage),
%         RT(i)               = time(firstPassage);
%         
%         switch boundPassed
%             case 1
%                 choice(i)           = 1;
%             case 2
%                 choice(i)           = -1;
%             otherwise
%                 disp(boundPassed);
%         end
%     end
% end

assert(length(unique(RT)) > 1);
RT = RT(~isnan(RT));
choice = choice(~isnan(choice));

% get the summary measure for each bin
qntls = quantile(RT, [0, 0.1, 0.3, 0.5, 0.7, 0.9, 1]);
rtbin =  discretize(RT, qntls);
pbias = splitapply(@nanmean, choice, rtbin);
plot([0, 0.1, 0.3, 0.5, 0.7, 0.9], pbias, 'color', color);


end