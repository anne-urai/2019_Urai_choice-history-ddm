

% ==================================================================
% now append all the eyelink files together
% ==================================================================

if ~isdeployed,
    addpath(genpath('~/code/RT_RDK'));
    addpath(genpath('~/code/Tools'));
    addpath('~/Documents/fieldtrip/');
    ft_defaults;
end

datapath = '~/Data/RT_RDK';
if ~exist('subjects', 'var'), subjects = [3:15 17:25]; end % 16 was excluded


for sj = subjects,
    clearvars -except subjects sj datapath;

    if exist(sprintf('%s/P%02d_alleye.mat', datapath, sj), 'file'),
        continue;
    end

    % check if the full dataset is not there yet
    cd(sprintf('%s/P%02d/pupil', datapath, sj));
    eyelinkfiles = dir(sprintf('P%02d*_eyeclean.mat', sj));

    % make sure these are in the right order!
    % otherwise, indexing of trials will go awry
    for f = 1:length(eyelinkfiles),
        scandat         = sscanf(eyelinkfiles(f).name, 'P%*d_s%d_b%d_*.mat');
        snum(f,:)       = scandat';
    end
    [sorted, sortidx]   = sort(snum(:,1)); % sort by session
    sorted(:,2)         = snum(sortidx, 2); % sort by block
    eyelinkfiles        = eyelinkfiles(sortidx);

    cfg             = [];
    cfg.inputfile   = {eyelinkfiles.name};
    cfg.outputfile  = sprintf('%s/P%02d_alleye.mat', datapath, sj);
    ft_appenddata(cfg);
end

% to save disk space
% eyelinkfiles = rdir(sprintf('%s/P%02d/pupil/P%02d*_eyeclean.mat', datapath, sj, sj));

% ==================================================================
% WRITE BEHAVIOURAL DATA AND SINGLE-TRIAL PUPIL MEASURES
% ==================================================================

clearvars -except subjects datapath
for sj = subjects,

    clearvars -except sj subjects alldat datapath;
    load(sprintf('%s/P%02d_alleye.mat', datapath, sj)); % contains all data over sessions
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

    assert(all(data.trialinfo(:,2)>50), 'wrongly coded stimulus');
    assert(all(data.trialinfo(:,4)>50), 'wrongly coded stimulus');

    % compute RT again
    RT          = data.trialinfo(:,6);

    % trialinfo matrix as it is
    baseline_pupil      = nan(length(data.trial), 1);
    response_pupil      = nan(length(data.trial), 1);
    decision_pupil      = nan(length(data.trial), 1);
    feedback_pupil      = nan(length(data.trial), 1);

    for t = 1:length(data.trial),

        % baseline
        stimstart = data.trialinfo(t,3)-data.trialinfo(t,1);
        baseline_pupil(t)  = mean(data.trial{t}(pupilchan, stimstart-0.5*data.fsample : stimstart));

        % dilation 1.5s after decision - subtract pretrial baseline
        responseonset  = data.trialinfo(t,7) - data.trialinfo(t,1);
        response_pupil(t) = mean(data.trial{t}(pupilchan, responseonset : responseonset+1.5*data.fsample)) - baseline_pupil(t);

        % dilation before feedback - subtract pretrial baseline
        fbonset  = data.trialinfo(t,9) - data.trialinfo(t,1);
        decision_pupil(t) = mean(data.trial{t}(pupilchan, fbonset - 0.5*data.fsample : fbonset)) - baseline_pupil(t);

        % dilation after feedback - subtract pretrial baseline
        feedback_pupil(t) = mean(data.trial{t}(pupilchan, fbonset : fbonset+1.5*data.fsample)) - baseline_pupil(t);

    end

    alldat{sj} = [data.trialinfo(:,2) data.trialinfo(:,4) RT data.trialinfo(:,5) ...
        data.trialinfo(:, [10 11 12]) sj*ones(size(RT)) ...
        baseline_pupil response_pupil decision_pupil feedback_pupil];

    switch sj
        case {3, 5, 15, 18}
            alldat{sj}(:, 7) = alldat{sj}(:, 7) - 1; % start at session 1
        case 17
            alldat{sj}(find(alldat{sj}(:, 7)) == 5, 7) = ...
                alldat{sj}(find(alldat{sj}(:, 7)) == 5, 7) - 1;
        case 20
            alldat{sj}(:, 7) = alldat{sj}(:, 7) - 2; % start at session 1
        case 4
            alldat{sj}(:, 7) = alldat{sj}(:, 7) - 2; % start at session 1
            assert(isequal(unique(alldat{sj}(:, 7))', 1:5), 'did not find 5 sessions');
        otherwise
            assert(isequal(unique(alldat{sj}(:, 7))', 1:5), 'did not find 5 sessions');
    end

    fprintf('\n subject %d, sessions %d %d %d %d %d %d %d \n\n', sj, unique(alldat{sj}(:, 7))');

    subplot(5,5,find(sj==subjects));
    plot(baseline_pupil, response_pupil, '.'); %l = lsline; set(l, 'color', 'k');
    axis tight; axis square; set(gca, 'tickdir', 'out');
end

disp('writing to table');
alldat2 = cat(1, alldat{:});

% write to csv for all subjects
t = array2table(alldat2, 'VariableNames', ...
    {'stimulus', 'response', 'rt', 'correct', ...
    'trialnr', 'blocknr', 'session', 'subj_idx',  ...
    'baseline_pupil', 'response_pupil', 'decision_pupil', 'feedback_pupil'});

% recode for HDDM
t.stimulus(t.stimulus == 90) = 0;
t.stimulus(t.stimulus == 270) = 1;
t.response(t.response == 90) = 0;
t.response(t.response == 270) = 1;

responseSigned = t.response;
responseSigned(responseSigned == 0) = -1;
t.prevresp  = circshift(responseSigned, 1);
t.prevrt    = circshift(log(t.rt), 1);
t.prevpupil = circshift(t.decision_pupil, 1);

% normalise pupil and rt within each block
nanzscore = @(x) (x - nanmean(x)) ./ nanstd(x);

for sj = unique(t.subj_idx)',
  for session = unique(t.session)',
    t.prevrt(find(t.subj_idx == sj & t.session == session)) = ...
    nanzscore(t.prevrt(find(t.subj_idx == sj & t.session == session)));
    t.prevpupil(find(t.subj_idx == sj & t.session == session)) = ...
    nanzscore(t.prevpupil(find(t.subj_idx == sj & t.session == session)));
  end
end

% remove trials where the previous trial was not immediately preceding
wrongtrls               = find([NaN; diff(t.trialnr)] ~= 1);
t.prevresp(wrongtrls)   = NaN;
t.prevrt(wrongtrls)     = NaN;
t.prevpupil(wrongtrls)  = NaN;

% correct some session nrss
t.session(t.session == 0) = 1;
t.session(t.subj_idx == 17 & t.session == 5) = 4;

writetable(t, sprintf('%s/HDDM/rtrdk_data_allsj.csv', datapath));
