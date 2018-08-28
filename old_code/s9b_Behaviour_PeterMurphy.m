function [rho] = s9b_Behaviour_PeterMurphy(subjects)
% replicate the patterns between baseline pupil and drift rate variability.
% prediction 1: higher bl pupil -> lower accuracy
% prediction 2: higher bl pupil -> greater difference between error and
% correct RTs
% from Murphy et al. 2014 PLoS Computational Biology
%
% Anne Urai, 2015

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
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT BASELINE VS DILATION
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    offset          = data_blcorr.trialinfo(:,6) - data_blcorr.trialinfo(:,1);
    prestim         = 0;
    poststim        = 1.5;
    decisiondilation  = shiftoffset_timelock(data_blcorr, [], offset, prestim, poststim, data_blcorr.fsample, 0);
    
    figure(fig1);
    subplot(5,5,find(sj==subjects));
    plot(bl, mean(decisiondilation.trial, 3), '.');
    axis tight; lsline;
    axis square;
    set(gca, 'box', 'on', 'tickdir', 'out');
    [r, pval] = corr(bl, mean(decisiondilation.trial, 3), 'type', 'Pearson');
    title(sprintf('rho %.2f p %.3f', r, pval));
    % waitforbuttonpress;
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FOR EACH SJ AND SESSION, COMPUTE THE ACCURACY BY 10 PUPIL BINS
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [binnedpupil(find(sj==subjects), :), binnedacc(find(sj==subjects), :), stdpup, stdacc] = ...
        divideintobins(bl, data.trialinfo(:,5), 10);
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FOR EACH SJ AND SESSION, COMPUTE THE ACCURACY BY 10 PUPIL BINS
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % use quantile rather than histcounts, we want each bin to contain the same nr of points!
    qs                                     = quantile(bl, 9);
    binnedRTcorr(find(sj==subjects), :)    = nan(1, 10);
    binnedRTerr(find(sj==subjects), :)     = nan(1, 10);
    
    for q = 1:length(qs) + 1,
        
        % determine which trials belong to this quantile
        if q == 1,
            findtrls = find(bl < qs(q));
        elseif q == length(qs) + 1,
            findtrls = find(bl >= qs(q-1));
        else
            findtrls = find(bl < qs(q) & bl >= qs(q-1));
        end
        
        % compute the mean pupil baseline in this bin
        binnedRTcorr(find(sj==subjects), q) = mean(RT(intersect(findtrls, find(data.trialinfo(:,5)==1))));
        binnedRTerr(find(sj==subjects), q) = mean(RT(intersect(findtrls, find(data.trialinfo(:,5)==0))));

    end
end

figure;
errorbar(1:10, mean(binnedpupil, std(binnedpupil))); xlabel('Pupil baseline'); ylabel('Accuracy');

end
