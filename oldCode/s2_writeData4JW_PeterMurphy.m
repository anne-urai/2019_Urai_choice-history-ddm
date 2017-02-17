
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WRITE BEHAVIOURAL DATA AND SINGLE-TRIAL PUPIL MEASURES
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;
cd ~/Data/PeterMurphy/

if ~exist('subjects', 'var'),
    % get all those subjects who have a folder
    cd ~/Data/PeterMurphy/;
    s = dir('P*');
    s = {s([s(:).isdir]).name};
    for i = 1:length(s), subjects(i) = str2num(s{i}(2:4)); end
end
ft_defaults;

for sj = subjects,
    
    clearvars -except sj subjects alldat;
    load(sprintf('~/Data/PeterMurphy/P%02d_alleye.mat', sj)); % contains all data over sessions
    data.fsample    = 250; % these data have not been resampled
    RT              = data.trialinfo(:,4);

    % remove RT outliers
    cfg                         = [];
    cfg.trials                  = true(1, length(data.trial));
    cfg.trials(find(RT > 2.99)) = false;
    cfg.trials(find(RT < 0.15)) = false;
    data            = ft_selectdata(cfg, data);
    pupilchan       = find(strcmp(data.label, 'EyePupil')==1);
    
    % compute RT again
    RT              = data.trialinfo(:,4);
    
    % trialinfo matrix as it is
    blpup       = nan(length(data.trial), 1);
    dilationpup = nan(length(data.trial), 1);
    
    for t = 1:length(data.trial),

        % baseline
        stimstart = data.trialinfo(t,3)-data.trialinfo(t,1);
        blpup(t)  = mean(data.trial{t}(find(data.time{t} < 0)));
       
    end
    
    alldat{sj} = [data.trialinfo(:,1) RT data.trialinfo(:,3) ...
        sj*ones(size(RT)) blpup ];
   
end

disp('writing to table');
alldat2 = cat(1, alldat{:});

% write to csv for all subjects
t = array2table(alldat2, 'VariableNames', ...
    {'stimulus',  'rt', 'response', ...
     'subj_idx', 'pupil_b'});
writetable(t, sprintf('petermurphy_data_allsj.csv'));
