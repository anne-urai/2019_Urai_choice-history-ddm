
ft_defaults;
if ~exist('subjects', 'var'), subjects = [3:15 17:25]; end
close all;
clc; clearvars -except subjects

set(0, 'DefaultAxesFontSize', 7);
fig1 = figure; %fig2 = figure; fig3 = figure;

for sj = subjects,
    
    clearvars -except sj subjects fig1 fig2 fig3 binnedpupil binnedacc binnedRTcorr binnedRTerr;
    clc;
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PREPARE
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    load(sprintf('~/Data/RT_DDM/P%02d_alleye.mat', sj)); % contains all data_blcorr over sessions
    data.fsample    = 100;
    pupilchan       = find(strcmp(data.label, 'EyePupil')==1);
    
    % remove noresp trials
    cfg        = [];
    cfg.trials = true(1, length(data.trial));
    cfg.trials(find(isnan(data.trialinfo(:,5)))) = false;
    data        = ft_selectdata(cfg, data);
    
    fprintf('removing %d out of %d trials \n', length(find(cfg.trials==0)), length(cfg.trials));
    RT = (data.trialinfo(:,6) - data.trialinfo(:,3)) ./ data.fsample;
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PRE-STIMULUS BASELINE APPLIED TO ALL DATA
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    disp('baseline correcting...');
    bl = nan(length(data.trial), 1);
    data_blcorr = data;
    
    for t = 1:length(data_blcorr.trial),
        
        startofref = data_blcorr.trialinfo(t,3) - data_blcorr.trialinfo(t,1);
        % use the 500ms before this (50 samples, resampled at 100Hz)
        try
            bl(t) = mean(data_blcorr.trial{t}(pupilchan, startofref-0.5*data_blcorr.fsample:startofref));
        catch
            % in case there are not enough samples before that
            bl(t) = mean(data_blcorr.trial{t}(pupilchan, 1:startofref));
        end
        data_blcorr.trial{t}(pupilchan, :) = data_blcorr.trial{t}(pupilchan, :) - bl(t); % subtract from whole trial
    end
   
    % find autocorrelation in the bl
    max_lag = 50;
    autocorr = xcorr(bl, 50);
    
end
