% SfN: perceptual learning

clear; close all; clc;
addpath(genpath('~/code/Tools'));
warning off;

set(groot, 'defaultaxesfontsize', 6, 'defaultaxestitlefontsizemultiplier', 1, ...
    'defaultaxestitlefontweight', 'bold', ...
    'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters');

% ============================================ %
% BEHAVIOURAL MEASURES OF LEARNING
% ============================================ %

usepath     = '~/Data/projects/0/neurodec/Data/MEG-PL';
results     = readtable(sprintf('%s/HDDM/summary/allindividualresults.csv', usepath));
subjects    = results.subjnr(results.session == 2);

% learning rate
t.driftrate_learn = results.v_stimulus_session__regress_dc_z_prevresp(results.session == 0 ...
    & ismember(results.subjnr, subjects));

% also compute dprime ratio

t.dprime_s1 = results.dprime(results.session == 1 & ismember(results.subjnr, subjects));
t.dprime_s2 = results.dprime(results.session == 2 & ismember(results.subjnr, subjects));
t.dprime_learn = t.dprime_s2 - t.dprime_s1;

% TO DO: get error bars from the posterior traces?

% ============================================ %
% GET MEG STUFF
% ============================================ %

% email Tobi 07.04.2017
% 1. comparison in the quality of sensory evidence encoding in visual cortex
% between session 1 and 2. Can be quantified with ROC analysis high vs low response distributions
% 2. comparison of the slope of mean beta-lateralization in motor cortex between S1 and 2

% ============================================ %
% 1. SENSORY EVIDENCE ENCODING
% ============================================ %

sensoryROC = nan(65, 2);

for s = 1:2,
    % get gamma-band response timecourses for every single trial
    load(sprintf('%s/GrandAverage/TOI/%s_S%d_dB_%sbaseline.mat', ...
        usepath, 'stim', s, 'common'));
    
    subjnr      = round(data.trialinfo(:, end) ./ 1000000);
    subjects    = unique(subjnr)'; subjects(isnan(subjects)) = [];
    for sj = subjects,
        
        thisidx     = find(sj == subjnr);
        thisdat     = data.avg(thisidx, :);
        thisinfo    = data.trialinfo(thisidx, :);
        
        % ROC BETWEEN STRONG AND WEAK TRIALS
        ref_timeidx  = 14:26;
        stim_timeidx = 63:75; % 150ms to 750ms from stimulus onset
        
        % the difference between stim and ref is what matters for the choice
        gammaBand_ref    = nanmean(thisdat(:, ref_timeidx), 2);
        gammaBand_stim   = nanmean(thisdat(:, stim_timeidx), 2);
        gammaBand        = gammaBand_stim - gammaBand_ref;
        
        % based on stimulus or response?
        rocval = rocAnalysis(gammaBand(thisinfo(:, 4) == -1), ...
            gammaBand(thisinfo(:, 4) == 1), 0, 0);
        
        sensoryROC(sj, s) = rocval.i;
    end
end

% subtract 0.5 for plotting
sensoryROC(isnan(mean(sensoryROC, 2)), :) = [];
sensoryROC = sensoryROC - 0.5;

t.visualgamma_s1 = sensoryROC(:, 1);
t.visualgamma_s2 = sensoryROC(:, 2);

% how much did this gamma-band ROC change from S1 to S2?
t.visualgamma_learn       = t.visualgamma_s2 - t.visualgamma_s1;

% ============================================ %
% 2. MOTOR BETA BUILD-UP
% ============================================ %

betabuildup = nan(65, 2);

for s = 1:2,
    % get gamma-band response timecourses for every single trial
    load(sprintf('%s/GrandAverage/TOI/%s_S%d_dB_%sbaseline.mat', ...
        usepath, 'resp', s, 'common'));
    
    subjnr      = round(data.trialinfo(:, end) ./ 1000000);
    subjects    = unique(subjnr)'; subjects(isnan(subjects)) = [];
    for sj = subjects,
        
        thisidx     = find(sj == subjnr);
        thisdat     = data.avg(thisidx, :);
        thisinfo    = data.trialinfo(thisidx, :);
        
        % compute lateralisation index by taking lefthand - righthand!
        beta    = nanmean(thisdat((thisinfo(:, 6) == 18), :), 1) - ...
            nanmean(thisdat((thisinfo(:, 6) == 12), :), 1);
        
        % response-locked beta-band build-up
        resp_timeidx = 135:149; % 150ms to 750ms from stimulus onset
        
        resp_timeidx = 63:75; % 150ms to 750ms from stimulus onset
        
        % FIT A SLOPE TO ALL TRIALS?
        % slope = polyfit(data.timename(resp_timeidx), beta(resp_timeidx), 1);
        %  bls     = regress(beta(resp_timeidx)', [ones(size(resp_timeidx))' data.timename(resp_timeidx)']);
        
        % robust regression against outliers
        brob    = robustfit(data.timename(resp_timeidx), beta(resp_timeidx));
        
        betabuildup(sj, s) = brob(2);
    end
end
betabuildup(isnan(mean(betabuildup, 2)), :) = [];

t.motorbeta_s1 = betabuildup(:, 1);
t.motorbeta_s2 = betabuildup(:, 2);

% how much did this gamma-band ROC change from S1 to S2?
t.motorbeta_learn       = t.motorbeta_s2 - t.motorbeta_s1;
figure;
corrplot(t, {'dprime_learn', 'driftrate_learn', 'visualgamma_learn', 'motorbeta_learn'});
print(gcf, '-dpdf', sprintf('%s/Figures/learning_improvements.pdf', usepath));

figure;
corrplot(t, {'dprime_s1', 'dprime_s2', 'visualgamma_s1', 'visualgamma_s2', 'motorbeta_s1', 'motorbeta_s2'});
print(gcf, '-dpdf', sprintf('%s/Figures/learning_allmeasures.pdf', usepath));


