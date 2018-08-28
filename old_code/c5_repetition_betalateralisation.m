
clear; close all; clc;
addpath(genpath('~/code/Tools'));
warning off;

usepath = '~/Data/projects/0/neurodec/Data/MEG-PL';

% behavioural results
results = readtable(sprintf('%s/HDDM/summary/allindividualresults.csv', usepath));

% compute repetition parameters from separate HDDM models
results.dc_rep__stimcoding_prevresp_dc = ...
    results.dc_1__stimcoding_prevresp_dc - results.dc_2__stimcoding_prevresp_dc;

results.z_rep__stimcoding_prevresp_z = ...
    results.z_2__stimcoding_prevresp_z - results.z_1__stimcoding_prevresp_z;

results.dc_rep__stimcoding_prevresp_dc_z = ...
    results.dc_1__stimcoding_prevresp_dc_z - results.dc_2__stimcoding_prevresp_dc_z;

results.z_rep__stimcoding_prevresp_dc_z = ...
    results.z_2__stimcoding_prevresp_dc_z - results.z_1__stimcoding_prevresp_dc_z;

params = {'criterionshift', ...
    'dc_rep__stimcoding_prevresp_dc', 'dc_rep__stimcoding_prevresp_dc_z',...
    'v_prevresp__regress_dc_prevresp','v_prevresp__regress_dc_z_prevresp', ...
    'z_rep__stimcoding_prevresp_z', 'z_rep__stimcoding_prevresp_dc_z', ...
    'z_prevresp__regress_z_prevresp', 'z_prevresp__regress_dc_z_prevresp'};

% ============================================ %
% ADD MEG PARAMETERS
% ============================================ %

cd('~/Dropbox/code/MEG');

% get table data, for next trial stuff
t           = readtable(sprintf('%s/CSV/2ifc_megdata_allsj.csv', usepath));
t.subjnr    = floor(t.idx/1e+006);
%t.repeat =
t.repeat = [NaN; (diff(t.resp) == 0)];
wrongTrls   = ([NaN; diff(t.trial)] ~= 1);
t.repeat(wrongTrls) = NaN;

what2plot.chans      = 'resp';
what2plot.contrast   = {'hand', 'repeat'}; % strong vs weak resp, has been flipped
what2plot.difference = {'hand'}; % strong vs weak resp, has been flipped

clf;
for s = 1:2,
    load(sprintf('%s/GrandAverage/TOI/%s_S%d_dB_%sbaseline.mat', ...
        usepath, 'resp', s, 'common'));
    [~, ~, ~, activity] = plotData(data, t, what2plot);
    close all;
    
    lateralisation = activity.avg(activity.contrast == 2, :) - activity.avg(activity.contrast == 1, :);
    lateralisation = nanmean(lateralisation(:, 11:26), 2); % samples 11 to 26 are during the reference interval
    
    % show
    hold on;
    plot(lateralisation, results.repetitioncorr(results.session == s), 'o');
    xlabel('Beta rebound previous trial');
    ylabel('P(repeat)'); lsline;
    hline(0, 'color', [0.5 0.5 0.5]);
    vline(0, 'color', [0.5 0.5 0.5]);
    plot(lateralisation, results.repetitioncorr(results.session == s), 'o');
    print(gcf, '-dpdf', sprintf('%s/HDDM/summary/repetition_corrplot_betaband_s%d.pdf', usepath, s));
    
end


% ============================================ %
% CORRELATION PLOT
% ============================================ %

% skip if some are missing
params = intersect(params, results.Properties.VariableNames');

% see how they correlate
corrplot(results, params);
print(gcf, '-dpdf', sprintf('~/Data/%s/HDDM/summary/repetition_corrplot_betaband.pdf', datasets{d}));

