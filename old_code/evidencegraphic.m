
close all;
cfg.timestep = 0.008;
cfg.time     = cfg.timestep:cfg.timestep:tmax;
cfg.a        = 1; % bound with, z = 0
cfg.cdW      = 0.1; % variance of normally distributed noise
cfg.v        = 2.5*cfg.timestep; % drift rate per timestep
cfg.z        = 0; % drift rate per timestep
%cfg.seed     = seed;
defcfg = cfg;

% generate a sample of evidence for each moment in time
evidence = cfg.v * ones(size(cfg.time));

% add an additional increase in the middle
evidence(round(length(evidence)*0.4):round(length(evidence)*0.7)) = ...
    evidence(round(length(evidence)*0.4):round(length(evidence)*0.7)) * 8;

% add noise
%rng(cfg.seed); % make sure this is the same each time
noise = randn(size(evidence)) .* cfg.cdW;
noise(1) = cfg.z; % start with no noise
noisyevidence = evidence + noise;

subplot(441);
gr = cbrewer('seq', 'Greens', 1);
plot(noisyevidence, 'color', gr(end, :) ); hold on;
plot(evidence, 'k');
axis tight; axis off;

subplot(442);
plot(cumsum(noisyevidence), 'color', gr(end, :));
axis tight; axis off;

integratedevidence = cumsum(noisyevidence);
y = integratedevidence;
firstPassage = find(y > cfg.a);
y(firstPassage:end) = NaN;
y(y > cfg.a) = NaN;
y(y < -cfg.a) = NaN;

print(gcf, '-depsc', '~/Dropbox/DonnerLab/graphics/integratedevidence.eps');
