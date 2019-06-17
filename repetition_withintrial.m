function repetition_withintrial()

% Code to fit the history-dependent drift diffusion models as described in
% Urai AE, de Gee JW, Tsetsos K, Donner TH (2019) Choice history biases subsequent evidence accumulation. eLife, in press.
%
% MIT License
% Copyright (c) Anne Urai, 2019
% anne.urai@gmail.com

% ========================================== %
% conditional response functions from White & Poldrack
% run on simulated rather than real data
% ========================================== %

addpath(genpath('~/code/Tools'));
warning off; % close all;
global datasets datasetnames mypath colors
groups = {'alternators', 'repeaters', 'all'};
groups = {'all'};

% ========================================== %
% START
% ========================================== %


for g = 1:length(groups),
    
    alldat.repeat = nan(6,4);
    alldat.bias = nan(6,4);
    alldat.rt = nan(6,4);
    
    
    close all; subplot(4,4,1); hold on;
    colors2 = cbrewer('qual', 'Set2', length(datasets));
    markers = {'d', 's', '^', 'v',  '>', '<'};
    alldata = {};

    for d = 1:length(datasets),
        
        filename    = dir(sprintf('%s/%s/*.csv', mypath, datasets{d}));
        data        = readtable(sprintf('%s/%s/%s', mypath, datasets{d}, filename.name));
        
        if contains(datasetnames{d}{2}, 'FD'),
            if d == 2,
                data.rt = data.rt - 0.25 + 0.75;
            elseif d == 3,
                data.rt = data.rt - 0.25 + 0.5;
            elseif d == 4,
                data.rt = data.rt - 0.25 + 0.75;
                
            end
        end
        
        % for those datasets with varying coherence, take only the difficult trials
        % FOR THE TWO DATASETS WITH VARYING COHERENCE, REMOVE THE HIGH LEVELS
        if any(ismember('coherence', data.Properties.VariableNames)),
            if max(data.coherence) == 81,
                data = data(data.coherence < 27, :);
            elseif max(data.coherence) == 0.3,
                data = data(data.coherence < 0.1, :);
            else % do nothing
            end
        end
        
        % divide RT into quantiles for each subject
        discretizeRTs = @(x) {discretize(x, [0 0.4 0.8 1.6 5])};
        discretizeRTs = @(x) {discretize(x, [0 quantile(x, [.2, .4, .6, .8, .95])])};

        % discretize into bins of RT
        rtbins = splitapply(discretizeRTs, data.rt, findgroups(data.subj_idx));
        data.rtbins = cat(1, rtbins{:});
        
        % ALSO COMPUTE PBIAS
        data.repeat = (sign(data.prevresp) == sign(data.response - 0.1));
        
        data.biased = data.repeat;
        [gr2, sjs] = findgroups(data.subj_idx);
        sjrep = splitapply(@nanmean, data.repeat, gr2);
        alternators = sjs(sjrep < 0.5);
        
        switch groups{g}
            case 'alternators'
                data(~ismember(data.subj_idx, alternators), :) = [];
            case 'repeaters'
                data(ismember(data.subj_idx, alternators), :) = [];
            otherwise
                % recode into biased and unbiased choices
                altIdx = ismember(data.subj_idx, alternators);
                data.biased(altIdx) = double(~(data.biased(altIdx))); % flip
        end
        
        if size(data, 1) == 0,
            continue;
        end
        
        % SPLIT REPETITION BIAS BY RT QUANTILES
        [gr, sjidx, rtbins]      = findgroups(data.subj_idx, data.rtbins);
        repetition               = array2table([sjidx, rtbins], 'variablenames', {'subj_idx', 'rtbin'});
        repetition.choice        = splitapply(@nanmean, data.repeat, gr); % choice proportion
        
        % make into a subjects by rtbin matrix
        mat_tmp = unstack(repetition, 'choice', 'rtbin');
        mat     = mat_tmp{:, 2:end}; % remove the last one, only has some weird tail
        
        bias               = array2table([sjidx, rtbins], 'variablenames', {'subj_idx', 'rtbin'});
        bias.choice        = splitapply(@nanmean, data.biased, gr); % choice proportion
        % make into a subjects by rtbin matrix
        mat_tmp           = unstack(bias, 'choice', 'rtbin');
        mat2              = mat_tmp{:, 2:end}; % remove the last one, only has some weird tail
        
        % also compute the mean RT for each subject and RT bin
        rtAvg               = array2table([sjidx, rtbins], 'variablenames', {'subj_idx', 'rtbin'});
        rtAvg.rt            = splitapply(@nanmean, data.rt, gr); % choice proportion
        xRTs                = unstack(rtAvg, 'rt', 'rtbin');
        xRTs                = xRTs{:, 2:end}; % remove the last one, only has some weird tail
        assert(isequal(size(mat), size(xRTs)), 'mismatch');
        
        % NOW PLOT

        errorbar(nanmean(xRTs), nanmean(mat), nanstd(mat) ./ sqrt(size(mat, 1)), '-', 'marker', markers{d},...
         'color', colors2(d, :), ...
            'capsize', 0, 'markerfacecolor', 'w',  'markersize', 3, 'markeredgecolor', colors2(d, :));

        % fixed effects across dataset
        data.subj_idx = data.subj_idx + 1000*d;
        alldata{end+1} = data(:, {'subj_idx', 'rt', 'repeat'});
    end

    end

    hline(0.5);

    data = cat(1, alldata{:});

    discretizeRTs = @(x) {discretize(x, [0 quantile(x, [0.1, .2, .4, .6, .8, .95])])};

    rtbins = splitapply(discretizeRTs, data.rt, findgroups(data.subj_idx));
    data.rtbins = cat(1, rtbins{:});
    
    % SPLIT REPETITION BIAS BY RT QUANTILES
    [gr, sjidx, rtbins]      = findgroups(data.subj_idx, data.rtbins);
    repetition               = array2table([sjidx, rtbins], 'variablenames', {'subj_idx', 'rtbin'});
    repetition.choice        = splitapply(@nanmean, data.repeat, gr); % choice proportion
    
    % make into a subjects by rtbin matrix
    mat_tmp = unstack(repetition, 'choice', 'rtbin');
    mat     = mat_tmp{:, 2:end}; % remove the last one, only has some weird tail
    
      % also compute the mean RT for each subject and RT bin
    rtAvg               = array2table([sjidx, rtbins], 'variablenames', {'subj_idx', 'rtbin'});
    rtAvg.rt            = splitapply(@nanmean, data.rt, gr); % choice proportion
    xRTs                = unstack(rtAvg, 'rt', 'rtbin');
    xRTs                = xRTs{:, 2:end}; % remove the last one, only has some weird tail
    assert(isequal(size(mat), size(xRTs)), 'mismatch');
    
    % ADD THE MEAN ACROSS DATASETS
    plot(nanmean(xRTs), nanmean(mat), '-', 'marker', '.', 'color', 'k');

    % set(gca, 'xtick', [0.4 0.8 1.6], 'xticklabelrotation', 45);
    % ylim([0.4 0.6]);
    axis tight;
    offsetAxes;
    ylabel('P(repeat)');
    %title(datasetnames{d});2
    %         if contains(datasetnames{d}{2}, 'RT'),
    xlabel('RT (s)');
    %         elseif contains(datasetnames{d}{2}, 'FD'),
    %             xlabel('RT from stim offset (s)');
    %         end
    set(gca, 'xcolor', 'k', 'ycolor', 'k');
    tightfig;
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/bias_withintrial_summary_quantiles.pdf', d)); % 3b
    disp('done repetition within trial')

    return;

    %%%%%%%%%%%%%%%%%
    % AVERAGE ACROSS DATASETs
    %%%%%%%%%%%%%%%%%
    
    close all; subplot(4,4,d); hold on;
    errorbar(nanmean(alldat.rt), nanmean(alldat.repeat), nanstd(alldat.repeat) ./ sqrt(size(alldat.repeat, 1)), '-bo', ...
        'capsize', 0, 'markerfacecolor', 'b', 'markeredgecolor', 'w');
    % errorbar(nanmean(alldat.rt), nanmean(alldat.bias), nanstd(alldat.bias) ./ sqrt(size(alldat.bias, 1)), '-ko', ...
    %     'capsize', 0, 'markerfacecolor', 'k', 'markeredgecolor', 'w');
    
    hline(0.5);
    set(gca, 'xtick', nanmean(alldat.rt), 'xticklabelrotation', 45);
    switch groups{g}
        case 'all'
           % ylim([0.48 0.56]);
    end
    
    offsetAxes;
    ylabel('P(repeat)');
    title(capitalize(groups{g}));
    xlabel('RT (s)');
    set(gca, 'xcolor', 'k', 'ycolor', 'k');
    tightfig;
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/repetition_withintrial_%s.pdf', groups{g})); % 3b
    
end

