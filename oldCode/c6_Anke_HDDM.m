% from the csv table, make an overview of repetition behaviour

% get a huge list with values for each participant
% can then work with this dataframe

clear; close all; clc;
addpath(genpath('~/code/Tools'));
warning off;

datasets = 'Anke_2afc_sequential/HDDM';
set(groot, 'defaultaxesfontsize', 8, 'defaultaxestitlefontsizemultiplier', 1, ...
    'defaultaxestitlefontweight', 'bold', ...
    'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters');
colormap viridis;

% ============================================ %
% INDIVIDUAL POINT ESTIMATES
% ============================================ %

clear; close;
dat = readtable('~/Data/Anke_2afc_sequential/HDDM/allindividualresults.csv');
% graphics
linspec = linspecer(3, 'qualitative');
% linspec = linspec([3 2 1], :); % blue first
subplot(331); hold on;
tp = unique(dat.transprob);
for t = 1:length(tp),
    plot(dat.v_prevresp__regress_dc_z_prevresp(dat.transprob == tp(t)), ...
        dat.v_prevpupil_prevresp__regress_dc_z_prevresp_prevpupil_prevrt(dat.transprob == tp(t)), ...
        '.', 'color', linspec(t, :));
end
xlabel('Repetition (dc)'); ylabel('Pupil x repetition (dc)'); 
axis square; axisNotSoTight; lsline;
for t = 1:length(tp),
    plot(dat.v_prevresp__regress_dc_z_prevresp(dat.transprob == tp(t)), ...
        dat.v_prevpupil_prevresp__regress_dc_z_prevresp_prevpupil_prevrt(dat.transprob == tp(t)), ...
        'o', 'markeredgecolor', 'w', 'markerfacecolor', linspec(t, :), 'markersize', 5);
end

subplot(332); hold on;
tp = unique(dat.transprob);
for t = 1:length(tp),
    plot(dat.v_prevresp__regress_dc_z_prevresp(dat.transprob == tp(t)), ...
        dat.v_prevrt_prevresp__regress_dc_z_prevresp_prevpupil_prevrt(dat.transprob == tp(t)), ...
        '.', 'color', linspec(t, :));
end
xlabel('Repetition (dc)'); ylabel('RT x repetition (dc)'); 
axis square; axisNotSoTight; lsline; 
for t = 1:length(tp),
    plot(dat.v_prevresp__regress_dc_z_prevresp(dat.transprob == tp(t)), ...
        dat.v_prevrt_prevresp__regress_dc_z_prevresp_prevpupil_prevrt(dat.transprob == tp(t)), ...
        'o', 'markeredgecolor', 'w', 'markerfacecolor', linspec(t, :), 'markersize', 5);
end

print(gcf, '-dpdf', sprintf('~/Data/%s/Figures/HDDMscatters.pdf', 'Anke_2afc_sequential/HDDM'));



% recode starting point values
z_link_func = @(x) 1 ./ (1 + exp(x));

% ============================================ %
% BETWEEN THE DIFFERENT TRANSITION PROBABILITY SESSIONS,
% HOW DO PEOPLE BEHAVE?
% ============================================ %

conditions = {'repetitive', 'neutral', 'alternating'};
vars = {'v_prevresp', 'v_prevrt_prevresp', 'v_prevpupil_prevresp', ...
    'z_prevresp', 'z_prevrt_prevresp', 'z_prevpupil_prevresp'};
close;
for v = 1:length(vars),
    subplot(3,3,v);
    hold on;
    for c = 1:length(conditions),
        dat = readtable(sprintf('~/Data/Anke_2afc_sequential/HDDM-%s/summary/regress_dc_z_prevresp_prevpupil_prevrt_all_traces_concat.csv', conditions{c}));
        if strcmp(vars{v}(1), 'z'),
            dat.(vars{v}) = z_link_func(dat.(vars{v}));
        end
        histogram(dat.(vars{v}), 'displaystyle', 'stairs', 'normalization', 'probability');
    end

    if strcmp(vars{v}(1), 'v'),
        vline(0, 'color', [0.5 0.5 0.5]);
    elseif strcmp(vars{v}(1), 'z'),
        vline(0.5, 'color', [0.5 0.5 0.5]);
    end
    xlabel(vars{v}, 'interpreter', 'none'); box off;
    
end

print(gcf, '-dpdf', sprintf('~/Data/%s/Figures/HDDMoverview.pdf', datasets));

