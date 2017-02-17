function [rho] = s9a_Pupil_RTvariability(subjects)
% from the individual pupil data

if ~exist('subjects', 'var'), subjects = [3:15 17:25]; end
close all;
clc; clearvars -except subjects
ft_defaults;

set(0, 'DefaultAxesFontSize', 7);
for sj = subjects,
    
    clearvars -except sj subjects rho binned
    clc; close all; figure;
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PREPARE
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    load(sprintf('~/Data/RT_DDM/P%02d_alleye.mat', sj)); % contains all data over sessions
    data.fsample    = 100;
    
    % remove noresp trials
    cfg             = [];
    cfg.trials      = true(1, length(data.trial));
    cfg.trials(find(isnan(data.trialinfo(:,5)))) = false;
    data            = ft_selectdata(cfg, data);
    pupilchan       = find(strcmp(data.label, 'EyePupil')==1);
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % RT VARIABILITY
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    subplot(4,4,1:2);
    RT          = (data.trialinfo(:,6)-data.trialinfo(:,3))/data.fsample;
    accuracy    = nan(size(RT));
    bin         = 50;
    nblocks     = 1:bin:length(data.trial);
    for b = nblocks,
        try
            accuracy(b) = nanmean(data.trialinfo(b:b+bin-1, 5));
        catch
            accuracy(b) = nanmean(data.trialinfo(b:end, 5));
        end
    end
    
    [hAx,hLine1,hLine2] = plotyy(1:length(data.trial), RT,  ...
        1:length(data.trial), accuracy);
    hAx(1).TickDir = 'out'; hAx(1).Box = 'off'; hAx(1).XLim = [0 length(data.trial)+1];
    hAx(2).TickDir = 'out'; hAx(2).Box = 'off'; hAx(2).XLim = [0 length(data.trial)+1];
    hLine1.LineStyle = 'none'; hLine1.Marker = '.';
    hLine2.LineStyle = 'none'; hLine2.Marker = '*';
    ylabel(hAx(1), 'Reaction time (s)');
    ylabel(hAx(2), 'Accuracy');
    xlabel('Trials');
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PRE-STIMULUS BASELINE APPLIED TO ALL DATA
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    disp('baseline correcting...');
    bl = nan(length(data.trial), 1);
    
    for t = 1:length(data.trial),
        
        startofref = data.trialinfo(t,3) - data.trialinfo(t,1);
        % use the 500ms before this (50 samples, resampled at 100Hz)
        try
            bl(t) = mean(data.trial{t}(pupilchan, startofref-0.5*data.fsample:startofref));
        catch
            % in case there are not enough samples before that
            bl(t) = mean(data.trial{t}(pupilchan, 1:startofref));
            assert(1==0)
        end
        data.trial{t}(pupilchan, :) = data.trial{t}(pupilchan, :) - bl(t); % subtract from whole trial
    end
    
    % add to trialinfo
    data.trialinfo(:, end+1) = bl;
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COMPUTE SINGLE TRIAL SCALARS FOR DECISION AND PUPIL
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    offset          = data.trialinfo(:,6) - data.trialinfo(:,1);
    prestim         = 0;
    poststim        = 1.5;
    decisionpupil   = shiftoffset_timelock(data, [], offset, prestim, poststim, data.fsample, 0);
    
    % decision-pupil, linear projection
    decpupmean = squeeze(mean(decisionpupil.trial));
    for t = 1:length(decisionpupil.trial),
       % data.trialinfo(t, 13) = linearprojection(squeeze(decisionpupil.trial(t, 1, :)), decpupmean);
       data.trialinfo(t, 13) = squeeze(mean(decisionpupil.trial(t, 1, :)));
    end
    
    offset          = data.trialinfo(:,8) - data.trialinfo(:,1);
    prestim         = 0;
    poststim        = 2;
    feedbackpupil   = shiftoffset_timelock(data, [], offset, prestim, poststim, data.fsample, 0);
    
    % feedback-pupil, linear projection
    fbpupmean = squeeze(mean(feedbackpupil.trial));
    for t = 1:length(feedbackpupil.trial),
        data.trialinfo(t, 14) = linearprojection(squeeze(feedbackpupil.trial(t, 1, :)), fbpupmean);
        data.trialinfo(t, 14) = squeeze(mean(feedbackpupil.trial(t, 1, :)));
    end
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT BASELINE VS DILATION
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    subplot(4,4,4);
    plot(data.trialinfo(:,12), data.trialinfo(:,13), '.');
    axis tight; lsline;
    axis square;
    set(gca, 'box', 'on', 'tickdir', 'out');
    [r, pval] = corr(bl, mean(decisionpupil.trial, 3), 'type', 'Pearson');
    xlabel({'pupil baseline'; '(% signal change)'});
    ylabel({'pupil dilation'; '(linear projection)'});
    title(sprintf('rho %.2f p %.3f', r, pval));
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % IN BINS OF TRIALS, COMPUTE THE RT VARIABILITY
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    bin     = 19; % how many trials per bin?
    
    RTvar   = nan(1, ceil(length(data.trial)/bin));
    decpup  = nan(1, ceil(length(data.trial)/bin));
    blpup   = nan(1, ceil(length(data.trial)/bin));
    fbpup   = nan(1, ceil(length(data.trial)/bin));
    acc     = nan(1, ceil(length(data.trial)/bin));
    ts      = 1:bin:length(data.trial);
    
    for t = ts,
        try
            % RT variability
            RTvar(find(t==ts))    = var(RT(t:t+bin-1));
        catch
            RTvar(find(t==ts))    = var(RT(t:end));
        end
        
        % pupil baseline mean
        try
            blpup(find(t==ts))   = mean(bl(t:t+bin-1));
        catch
            blpup(find(t==ts))   = mean(bl(t:end));
        end
        
        % pupil dilation mean
        try
            decpup(find(t==ts))   = mean(data.trialinfo(t:t+bin-1, 13));
        catch
            decpup(find(t==ts))   = mean(data.trialinfo(t:end, 13));
        end
        
        % feedback  mean
        try
            fbpup(find(t==ts))   = mean(data.trialinfo(t:t+bin-1, 14));
            % add the accuracy
            acc(find(t==ts))      = mean(data.trialinfo(t:t+bin-1, 5));
            
        catch
            fbpup(find(t==ts))   = mean(data.trialinfo(t:end, 14));
            acc(find(t==ts))     = mean(data.trialinfo(t:end, 5));
            
        end
    end
    assert(~any(isnan(RTvar)));
    assert(~any(isnan(decpup)));
    assert(~any(isnan(blpup)));
    
    subplot(4,4,9);
    plot(blpup, RTvar, '.'); axis tight; set(gca, 'tickdir', 'out');
    [rho.baseline(find(sj==subjects)), pval] = corr(RTvar', blpup', 'type', 'Spearman');
    title(sprintf('rho %.2f p %.3f', rho.baseline(find(sj==subjects)), pval));
    lsline; xlabel({'Baseline pupil'; '(% signal change)'}); ylabel('RT variance');
    
    subplot(4, 4, 10);
    plot(decpup, RTvar, '.'); axis tight; set(gca, 'tickdir', 'out');
    [rho.decision(find(sj==subjects)), pval] = corr(RTvar', decpup', 'type', 'Spearman');
    title(sprintf('rho %.2f p %.3f', rho.decision(find(sj==subjects)), pval));
    lsline; xlabel({'Decision pupil'; '(linear projection)'});
    
    subplot(4, 4, 11);
    plot(fbpup, RTvar, '.'); axis tight; set(gca, 'tickdir', 'out');
    [rho.feedback(find(sj==subjects)), pval] = corr(RTvar', fbpup', 'type', 'Spearman');
    title(sprintf('rho %.2f p %.3f', rho.feedback(find(sj==subjects)), pval));
    lsline; xlabel({'Feedback pupil'; '(linear projection)'});
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SAVE BINNED DATA ACROSS SUBJECTS
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    bin = 10; % divide everything in 10 bins
    
    % the 3rd and 4th output args from divideintobins are std
    [binned.baseline.pupil(find(sj==subjects), :), ~, ...
        ~, binned.baseline.RT(find(sj==subjects), :)] ...
        = divideintobins(data.trialinfo(:, 12), RT, bin);
    
    [binned.decision.pupil(find(sj==subjects), :), ~, ...
        ~, binned.decision.RT(find(sj==subjects), :)] ...
        = divideintobins(data.trialinfo(:,13), RT, bin);
    
    [binned.feedback.pupil(find(sj==subjects), :), ~, ...
        ~, binned.feedback.RT(find(sj==subjects), :)] ...
        = divideintobins(data.trialinfo(:,14), RT, bin);
    
    % as a control, compute the accuracy for each fb pupil bin and the RT
    [binned.feedback.pupil(find(sj==subjects), :), ...
        binned.feedback.accuracy(find(sj==subjects), :)] ...
        = divideintobins(data.trialinfo(:,14), data.trialinfo(:,5), bin);
    
    suptitle(sprintf('P%02d', sj));
    set(gcf, 'PaperPositionMode', 'auto');
    print(gcf, '-dpdf', '-painters', sprintf('~/Data/RT_DDM/Figures/P%02d_allpupilRT.pdf', sj));
    
end

save('~/Data/RT_DDM/allRTvariabibility.mat', 'rho', 'binned', 'subjects', 'bin');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT GRAND AVERAGE
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;
load('~/Data/RT_DDM/allRTvariabibility.mat');

figure;
subplot(3,3,1);
ploterr(mean(binned.baseline.pupil), mean(binned.baseline.RT), ...
    std(binned.baseline.pupil) / sqrt(length(subjects)), ...
    std(binned.baseline.RT) / sqrt(length(subjects)), 'k*', 'hhxy', 0.1);
set(gca, 'tickdir', 'out', 'box', 'off'); xlabel('Baseline pupil'); ylabel('RT variance');

subplot(3,3,2);
ploterr(mean(binned.decision.pupil), mean(binned.decision.RT), ...
    std(binned.decision.pupil) / sqrt(length(subjects)), ...
    std(binned.decision.RT) / sqrt(length(subjects)), 'k*', 'hhxy', 0.1);
set(gca, 'tickdir', 'out', 'box', 'off'); xlabel('Decision pupil'); %ylabel('RT variance');

hold on;
plot(binned.decision.pupil', binned.decision.RT', '.-');
axis tight;

subplot(3,3,3);
ploterr(mean(binned.feedback.pupil), mean(binned.feedback.RT), ...
    std(binned.feedback.pupil) / sqrt(length(subjects)), ...
    std(binned.feedback.RT) / sqrt(length(subjects)), 'k*', 'hhxy', 0.1);
set(gca, 'tickdir', 'out', 'box', 'off'); xlabel('Feedback pupil'); %ylabel('RT variance');

subplot(3,3,4);
ploterr(mean(binned.feedback.pupil), mean(binned.feedback.accuracy), ...
    std(binned.feedback.pupil) / sqrt(length(subjects)), ...
    std(binned.feedback.accuracy) / sqrt(length(subjects)), 'k*', 'hhxy', 0.1);
set(gca, 'tickdir', 'out', 'box', 'off'); xlabel('Feedback pupil'); ylabel('Accuracy');


suptitle(sprintf('Grand Average, n = %d, bins of %d trials', length(subjects), bin));
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, '-dpdf', '-painters', '~/Data/RT_DDM/Figures/GA_binned_PupilbyRTvar.pdf');


end
