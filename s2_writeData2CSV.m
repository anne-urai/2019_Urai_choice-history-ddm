
if ~isdeployed,
  addpath(genpath('~/code/RT_RDK'));
  addpath(genpath('~/code/Tools'));
  addpath('~/Documents/fieldtrip/');
  ft_defaults;
end

datapath = '~/Data/RT_RDK';
if ~exist('subjects', 'var'), subjects = [8:15 17:25]; end % 16 was excluded

for sj = subjects,
  clearvars -except subjects sj datapath;

    % ==================================================================
    % now append all the eyelink files together
    % ==================================================================

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

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CORRECT SESSION NRS FOR SOME PARTICIPANTS
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0,
    clearvars -except sj subjects
    load(sprintf('~/Data/RT_DDM/P%02d_alleye.mat', sj));

    sessionsdone = unique(data.trialinfo(:, end));
    switch sj
    case 4
      assert(isequal(sessionsdone', 1:7), 'did not find 5 sessions');
    case 17
      data.trialinfo(find(data.trialinfo(:, end)==5), end) = 4;
      sessionsdone = unique(data.trialinfo(:, end));
      assert(isequal(sessionsdone', 1:4), 'did not find 5 sessions');
    case 18
      data.trialinfo(:, end) = data.trialinfo(:, end) - 1;
      sessionsdone = unique(data.trialinfo(:, end));
      assert(isequal(sessionsdone', 1:4), 'did not find 5 sessions');
    case 20
      data.trialinfo(:, end) = data.trialinfo(:, end) - 1;
      sessionsdone = unique(data.trialinfo(:, end));
      assert(isequal(sessionsdone', 1:4), 'did not find 5 sessions');
    otherwise
      assert(isequal(sessionsdone', 1:5), 'did not find 5 sessions');
    end

    save(sprintf('%s/P%02d_alleye.mat', datapath, sj), 'data');
  end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WRITE BEHAVIOURAL DATA AND SINGLE-TRIAL PUPIL MEASURES
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
