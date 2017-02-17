
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WRITE BEHAVIOURAL DATA AND SINGLE-TRIAL PUPIL MEASURES
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;
cd ~/Data/RT_DDM/
if ~exist('subjects', 'var'), subjects = [3:15 17:25]; end % 16 was excluded
ft_defaults;

for sj = subjects,
    
    clearvars -except sj subjects alldat;
    load(sprintf('~/Data/RT_DDM/P%02d_alleye.mat', sj)); % contains all data over sessions
    data.fsample    = 100;
    RT              = data.trialinfo(:,6);
    
    % remove trials without a response
    cfg             = [];
    cfg.trials      = true(1, length(data.trial));
    cfg.trials(find(isnan(data.trialinfo(:,4)))) = false;
    
    % remove RT outliers
    cfg.trials(find(RT > 2.99)) = false;
    cfg.trials(find(RT < 0.15)) = false;
    data            = ft_selectdata(cfg, data);
    pupilchan       = find(strcmp(data.label, 'EyePupil')==1);
    
    % correct weirdness
    data.trialinfo(find(data.trialinfo(:,2)==27), 2) = 270;
    data.trialinfo(find(data.trialinfo(:,2)==9), 2) = 90;
    data.trialinfo(find(data.trialinfo(:,4)==27), 4) = 270;
    data.trialinfo(find(data.trialinfo(:,4)==9), 4) = 90;
    
    assert(all(data.trialinfo(:,2)>50), 'wrongly coded stimulus');
    assert(all(data.trialinfo(:,4)>50), 'wrongly coded stimulus');
    
    % compute RT again
    RT          = data.trialinfo(:,6);
    
    % trialinfo matrix as it is
    blpup       = nan(length(data.trial), 1);
    dilationpup = nan(length(data.trial), 1);
    fbpup       = nan(length(data.trial), 1);
    
    for t = 1:length(data.trial),

        % baseline
        stimstart = data.trialinfo(t,3)-data.trialinfo(t,1);
        blpup(t)  = mean(data.trial{t}(pupilchan, stimstart-0.5*data.fsample : stimstart));
        
        % dilation after decision - subtract pretrial baseline
        responseonset  = data.trialinfo(t,7) - data.trialinfo(t,1);
        dilationpup(t) = mean(data.trial{t}(pupilchan, responseonset : responseonset+1.5*data.fsample)) - blpup(t);
        
        % dilation after feedback - subtract pretrial baseline
        fbonset  = data.trialinfo(t,9) - data.trialinfo(t,1);
        fbpup(t) = mean(data.trial{t}(pupilchan, fbonset : fbonset+1.5*data.fsample)) - blpup(t);
        
    end
    
    alldat{sj} = [data.trialinfo(:,2) data.trialinfo(:,4) RT data.trialinfo(:,5) ...
        data.trialinfo(:, [10 11 12]) sj*ones(size(RT)) ...
        blpup dilationpup fbpup];
    
    fprintf('\n subject %d, sessions %d %d %d %d %d %d %d \n', sj, unique(alldat{sj}(:, 7))');
    
    subplot(5,5,find(sj==subjects));
    plot(blpup, dilationpup, '.'); l = lsline; set(l, 'color', 'k');
    axis tight; axis square; set(gca, 'tickdir', 'out');
end

disp('writing to table');
alldat2 = cat(1, alldat{:});

% write to csv for all subjects
t = array2table(alldat2, 'VariableNames', ...
    {'stimulus', 'choice',  'rt', 'response', ...
    'trialnr', 'blocknr', 'sessionnr', 'subj_idx',  ...
  'pupil_b', 'pupil', 'pupil_fb'});
writetable(t, sprintf('fleurs_data_allsj.csv'));
