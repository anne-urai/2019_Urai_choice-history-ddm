function results = b3b_behaviouralMetrics(alldata)

% Code to fit the history-dependent drift diffusion models described in
% Urai AE, Gee JW de, Donner TH (2018) Choice history biases subsequent evidence accumulation. bioRxiv:251595
%
% MIT License
% Copyright (c) Anne Urai, 2018
% anne.urai@gmail.com

warning off;

% ========================================== %
% COMPUTE FOR ALL SUBJECTS
% ========================================== %

% preallocate variables
varnames = {'subjnr', 'session', 'dprime', 'accuracy', 'criterion', 'bias', 'abscriterion', 'rt', ...
    'criterionshift', 'criterionshift_prct', ...
    'repetition_prevcorrect', 'repetition_preverror', ...
    'pupil_correct', 'pupil_error', 'rt_correct', 'rt_error', ...
    'rt_valid_slow_correct', 'rt_valid_fast_correct', 'rt_invalid_slow_correct', 'rt_invalid_fast_correct', ...
    'rt_valid_slow_error', 'rt_valid_fast_error', 'rt_invalid_slow_error', 'rt_invalid_fast_error', ...
    'accuracy_valid_slow', 'accuracy_valid_fast', 'accuracy_invalid_slow', 'accuracy_invalid_fast'};

nrSess          = length(unique(alldata.session)) + 1;
results         = array2table(nan(length(unique(alldata.subj_idx))*nrSess, length(varnames)), 'variablenames', varnames);
results.drug    = repmat({'NaN'}, length(unique(alldata.subj_idx))*nrSess, 1);
results.criterionshift_prct = nan(height(results), 2); % will contain two points for error bars


if any(~cellfun(@isempty, strfind(alldata.Properties.VariableNames, 'transitionprob'))),
    results.criterionshift_neutral = nan(size(results.criterionshift));
    results.criterionshift_alternating = nan(size(results.criterionshift));
    results.criterionshift_repetitive = nan(size(results.criterionshift));
    
    results.repetition_neutral = nan(size(results.criterionshift));
    results.repetition_alternating = nan(size(results.criterionshift));
    results.repetition_repetitive = nan(size(results.criterionshift));
end

% preallocate dprime for different coherence levels
if sum(strcmp(alldata.Properties.VariableNames, 'coherence')) > 0,
    cohlevels = unique(alldata.coherence);
    for c = 1:length(cohlevels),
        vrnm = ['dprime_c' num2str(cohlevels(c)*100)];
        vrnm = regexprep(vrnm, '\.', '\_'); % replace points in varname
        results.(vrnm) = nan(size(results.dprime));
    end
end

% measures that are modulated by previous trial RT or pupil
metrics = {'dprime', 'criterion', 'abscriterion', 'accuracy', 'repetition', ...
    'stimrepetition', 'repetitioncrit', 'criterionshift', 'handshift', ...
    'nextdprime', 'nextabscriterion'};
modulation = false;

if modulation,
    
    % modulation factors
    factors = {'pupil', 'rt', 'pupil_correct', 'pupil_error', 'rt_correct', 'rt_error'};
    for m = 1:length(metrics),
        % first, without modulation
        results.([metrics{m}]) = nan(size(results.subjnr));
        
        % then modulated by each of the factors
        for f = 1:length(factors),
            results.([factors{f} '_' metrics{m} '_corr']) = nan(size(results.subjnr));
            for b = 1:3, % binned
                results.(sprintf('%s_%s_bin%d', metrics{m}, factors{f}, b)) = nan(size(results.subjnr));
            end
        end
    end
end

% get all data
subjects      = unique(alldata.subj_idx)';
if modulation,
    alldata.pupil = alldata.decision_pupil;
end

if all(cellfun(@isempty, strfind(alldata.Properties.VariableNames, 'prevrt'))),
    % normalise previous RT per block
    blockchange = find(diff(alldata.block) > 0);
    blocknrs = zeros(height(alldata), 1);
    for b = 1:length(blockchange)-1,
        blocknrs(blockchange(b)+1:blockchange(b+1)) = blocknrs(blockchange(b))+1;
    end
    blocknrs(blockchange(end)+1:end) = blocknrs(blockchange(end))+1;
    for b = unique(blocknrs)',
        alldata.prevrt(blocknrs == b) = nanzscore(log(alldata.prevrt(blocknrs == b)));
    end
end

% recode correct
if all(cellfun(@isempty, strfind(alldata.Properties.VariableNames, 'correct'))),
    tmpstim = (alldata.stimulus > 0);
    tmpresp = (alldata.response > 0);
    alldata.correct = (tmpstim == tmpresp);
end

% only MEG-PL data has starthand
if isfield(alldata, 'startHand'),
    alldata.startHand(alldata.startHand > 20) = nan;
else
    alldata.startHand = nan(size(alldata.subj_idx));
end

% for criterion shift
alldata.nextstim = circshift(alldata.stimulus, -1);
alldata.nextresp = circshift(alldata.response, -1);
try
    alldata.nextstim((diff(alldata.trial) ~= 1)) = NaN;
    alldata.nextresp((diff(alldata.trial) ~= 1)) = NaN;
end

% for mulder et al. analysis
alldata.prevstim = circshift(alldata.stimulus, 1);
alldata.prevresp = circshift(alldata.response, 1);

try
    wrongtrls        = find([NaN; diff(alldata.trial)] ~= 1);
    alldata.prevstim(wrongtrls) = NaN;
    alldata.prevresp(wrongtrls) = NaN;
end

% ========================================== %
% READY, SET, GO
% ========================================== %

icnt = 0;
for sj = subjects,
    for s = [0 unique(alldata.session)'],
        
        icnt                    = icnt + 1;
        results.subjnr(icnt)    = sj;
        results.session(icnt)   = s;
        
        switch s
            case 0
                % all sessions together
                data                    = alldata(alldata.subj_idx == sj, :);
            otherwise % split by session
                data                    = alldata(alldata.subj_idx == sj & alldata.session == s, :);
        end
        data(isnan(data.response), :) = [];
        
        % some people don't have pupil data in each session
        if isempty(data),
            fprintf('skipping sj %d, session %d \n', sj, s);
            continue;
        end
        
        % ========================================== %
        % GENERAL STUFF
        % ========================================== %
        
        [d, c] = dprime(data.stimulus, data.response);
        results.dprime(icnt)        = d;
        assert(~isnan(d), 'dprime cannot be NaN');
        results.criterion(icnt)     = c;
        results.abscriterion(icnt)  = abs(c);
        results.accuracy(icnt)      = nanmean(data.correct);
        results.rt(icnt)            = nanmedian(data.rt);
        results.bias(icnt)          = nanmean(data.response);
        
        if sum(strcmp(data.Properties.VariableNames, 'coherence')) > 0,
            cohlevels = unique(data.coherence);
            for c = 1:length(cohlevels),
                vrnm = ['dprime_c' num2str(cohlevels(c)*100)];
                vrnm = regexprep(vrnm, '\.', '\_'); % replace points in varname
                % disp(vrnm);
                results.(vrnm)(icnt) = ...
                    dprime(data.stimulus(data.coherence == cohlevels(c)), ...
                    data.response(data.coherence == cohlevels(c)));
                assert(~isnan(results.(vrnm)(icnt)));
            end
        end
        
        % measure of repetition behaviour
        % data.repeat = [~(abs(diff(data.response)) > 0); NaN];
        % data.stimrepeat = [~(abs(diff(data.stimulus)) > 0); NaN];
        
        % 01.10.2017, use the same metric as in MEG, A1c_writeCSV.m
        data.repeat = [NaN; (diff(data.response) == 0)];
        data.stimrepeat = [NaN; (diff(data.response) == 0)];
        wrongTrls   = ([NaN; diff(data.trial)] ~= 1);
        data.repeat(wrongTrls) = NaN;
        data.stimrepeat(wrongTrls) = NaN;
        
        results.repetition(icnt)        = nanmean(data.repeat);
        results.stimrepetition(icnt)    = nanmean(data.stimrepeat);
        
        % also compute this after error and correct trials
        results.repetition_prevcorrect(icnt) = nanmean(data.repeat((data.prevstim > 0) == (data.prevresp > 0)));
        results.repetition_prevcorrect(icnt) = nanmean(data.repeat((data.prevstim > 0) ~= (data.prevresp > 0)));
        
        % criterion based on repetition and stimulus sequences
        [~, c] = dprime(data.stimrepeat, data.repeat);
        results.repetitioncrit(icnt)    = -c;
        
        % criterion based on next trial bias, then collapsed
        results.criterionshift(icnt)    = criterionshift(data.response, data.nextstim, data.nextresp);
        
        % add: a bootstrapped measure of criterionshift for error bars in correlation plot
        % 05.09.2017, after talk with JW - this takes forever!
        %bootstat = bootstrp(1000, @criterionshift,data.response,data.nextstim, data.nextresp);
        %results.criterionshift_prct(icnt, :) = prctile(bootstat, [0.25 0.75]);
        
        if s == 0,
            thispersonsbias = results.repetition(icnt) - results.stimrepetition(icnt);
        end
        
        % does the random hand they press cause a bias?
        results.handshift(icnt)         = criterionshift(data.startHand, data.stimulus, data.response);
        
        try
            % pupilstuff
            results.pupil_error(icnt)       = nanmean(data.pupil(data.correct == 0));
            results.pupil_correct(icnt)     = nanmean(data.pupil(data.correct == 1));
        end
        
        results.rt_error(icnt)          = nanmedian(data.rt(data.correct == 0));
        results.rt_correct(icnt)        = nanmedian(data.rt(data.correct == 1));
        
        % ========================================== %
        %   for data with transition probability, compute criterionshift separately
        % ========================================== %
        
        if s == 0 & any(~cellfun(@isempty, strfind(data.Properties.VariableNames, 'transitionprob'))),
            
            results.criterionshift_neutral(icnt) = criterionshift(data.response(data.transitionprob == 0.5), ...
                data.nextstim(data.transitionprob == 0.5), data.nextresp(data.transitionprob == 0.5));
            results.criterionshift_alternating(icnt) = criterionshift(data.response(data.transitionprob == 0.2), ...
                data.nextstim(data.transitionprob == 0.2), data.nextresp(data.transitionprob == 0.2));
            results.criterionshift_repetitive(icnt) = criterionshift(data.response(data.transitionprob == 0.8), ...
                data.nextstim(data.transitionprob == 0.8), data.nextresp(data.transitionprob == 0.8));
            
            results.repetition_neutral(icnt) = nanmean(data.repeat(data.transitionprob == 0.5));
            results.repetition_alternating(icnt) = nanmean(data.repeat(data.transitionprob == 0.2));
            results.repetition_repetitive(icnt) = nanmean(data.repeat(data.transitionprob == 0.8));
            
        end
        
        % ========================================== %
        % MULDER ET AL. 2012
        % http://www.jneurosci.org/content/32/7/2335.long
        % treat the previous choice as a cue
        % ========================================== %
        
        if thispersonsbias < 0, % alternators
            % only previous trials that are correct
            validtrls     = ((data.prevresp == data.prevstim) & (data.prevresp ~= data.stimulus));
            invalidtrls   = ((data.prevresp == data.prevstim) & (data.prevresp == data.stimulus));
            %neutraltrls   = ones(height(data), 1); % only previous correct trls
        elseif thispersonsbias > 0, % repeaters
            % only previous trials that are correct
            validtrls     = ((data.prevresp == data.prevstim) & (data.prevresp == data.stimulus));
            invalidtrls   = ((data.prevresp == data.prevstim) & (data.prevresp ~= data.stimulus));
            %neutraltrls   = ones(height(data), 1); % only previous correct trls
        else
            validtrls       = ones(height(data), 1);
            invalidtrls     = ones(height(data), 1);
            %neutraltrls    = ones(height(data), 1)
        end
        
        % data.rt = nanzscore(log(data.rt));
        results.rt_invalid_fast_correct(icnt)  = nanmedian(data.rt(invalidtrls & data.correct == 1 & ...
            data.prevrt < median(data.prevrt(invalidtrls & data.correct == 1))));
        results.rt_invalid_slow_correct(icnt)  = nanmedian(data.rt(invalidtrls & data.correct == 1 & ...
            data.prevrt > median(data.prevrt(invalidtrls & data.correct == 1))));
        results.rt_valid_fast_correct(icnt)    = nanmedian(data.rt(validtrls & data.correct == 1 & ...
            data.prevrt < median(data.prevrt(validtrls & data.correct == 1))));
        results.rt_valid_slow_correct(icnt)    = nanmedian(data.rt(validtrls & data.correct == 1 & ...
            data.prevrt > median(data.prevrt(validtrls & data.correct == 1))));
        
        results.rt_invalid_fast_error(icnt)  = nanmedian(data.rt(invalidtrls & data.correct == 0 & ...
            data.prevrt < median(data.prevrt(invalidtrls & data.correct == 0))));
        results.rt_invalid_slow_error(icnt)  = nanmedian(data.rt(invalidtrls & data.correct == 0 & ...
            data.prevrt > median(data.prevrt(invalidtrls & data.correct == 0))));
        results.rt_valid_fast_error(icnt)    = nanmedian(data.rt(validtrls & data.correct == 0 & ...
            data.prevrt < median(data.prevrt(validtrls & data.correct == 0))));
        results.rt_valid_slow_error(icnt)    = nanmedian(data.rt(validtrls & data.correct ==0 & ...
            data.prevrt > median(data.prevrt(validtrls & data.correct == 0))));
        
        results.accuracy_invalid_fast(icnt)    = nanmean(data.correct(invalidtrls & ...
            data.prevrt < median(data.prevrt(invalidtrls))));
        results.accuracy_invalid_slow(icnt)    = nanmean(data.correct(invalidtrls & ...
            data.prevrt > median(data.prevrt(invalidtrls))));
        results.accuracy_valid_fast(icnt)    = nanmean(data.correct(validtrls & ...
            data.prevrt < median(data.prevrt(validtrls))));
        results.accuracy_valid_slow(icnt)    = nanmean(data.correct(validtrls & ...
            data.prevrt > median(data.prevrt(validtrls))));
        
        % ========================================== %
        % DPRIME, CRITERION, ABSOLUTE CRITERION
        % correlation to pupil and RT
        % ========================================== %
        
        if modulation,
            flds = {'pupil', 'rt', 'pupil_correct', 'pupil_error', 'rt_correct', 'rt_error'};
            nbins = 8; % with 5 bins, correlation quite unstable
            for f = 1:length(flds),
                
                % for error and correct, make a temporary field
                if strfind(flds{f}, '_correct'),
                    useFld = data.(flds{f}(1:end-8));
                    useFld(data.correct == 0) = NaN;
                elseif strfind(flds{f}, '_error'),
                    useFld = data.(flds{f}(1:end-6));
                    useFld(data.correct == 1) = NaN;
                else
                    useFld = data.(flds{f});
                end
                
                % sanity check for those without pupil stuff
                if all(isnan(useFld)), continue; end
                
                % ========================================== %
                % correlate across five bins
                % ========================================== %
                
                binIdx      = discretize(useFld, [-inf quantile(useFld, nbins-1) inf]);
                try % will put NaNs if not enough trials
                    a           = splitapply(@nanmean, data.correct, binIdx);
                    [~, r]      = splitapply(@dprime, data.stimrepeat, data.repeat, binIdx);
                    rep         = splitapply(@nanmean, data.repeat, binIdx);
                    stimrep     = splitapply(@nanmean, data.stimrepeat, binIdx);
                    cs          = splitapply(@criterionshift, data.response, data.nextstim, data.nextresp, binIdx);
                    [nextd, nextc]      = splitapply(@dprime, data.nextstim, data.nextresp, binIdx);
                    
                    % see corrFunc defined below
                    results.([flds{f} '_dprime_corr'])(icnt)            = corrFunc(d);
                    results.([flds{f} '_criterion_corr'])(icnt)         = corrFunc(c);
                    results.([flds{f} '_abscriterion_corr'])(icnt)      = corrFunc(abs(c));
                    results.([flds{f} '_accuracy_corr'])(icnt)          = corrFunc(a);
                    results.([flds{f} '_repetition_corr'])(icnt)        = corrFunc(rep);
                    results.([flds{f} '_stimrepetition_corr'])(icnt)    = corrFunc(stimrep);
                    results.([flds{f} '_repetitioncrit_corr'])(icnt)    = corrFunc(-r);
                    results.([flds{f} '_criterionshift_corr'])(icnt)    = corrFunc(cs);
                    results.([flds{f} '_nextdprime_corr'])(icnt)        = corrFunc(nextd);
                    results.([flds{f} '_nextabscriterion_corr'])(icnt)  = corrFunc(abs(nextc));
                end
                
                % ========================================== %
                % divide into 3 bins for viz
                % ========================================== %
                
                binIdx = discretize(useFld, [-inf quantile(useFld, 2) inf]);
                
                [d, c] = splitapply(@dprime, data.stimulus, data.response, binIdx);
                a      = splitapply(@nanmean, data.correct, binIdx);
                [~, r] = splitapply(@dprime, data.stimrepeat, data.repeat, binIdx);
                rep    = splitapply(@nanmean, data.repeat, binIdx);
                cs     = splitapply(@criterionshift, data.response, data.nextstim, data.nextresp, binIdx);
                [nextd, nextc]      = splitapply(@dprime, data.nextstim, data.nextresp, binIdx);
                
                for b = 1:length(d),
                    
                    results.(sprintf('dprime_%s_bin%d', flds{f}, b))(icnt)      = d(b);
                    results.(sprintf('criterion_%s_bin%d', flds{f}, b))(icnt)   = c(b);
                    results.(sprintf('abscriterion_%s_bin%d', flds{f}, b))(icnt)   = abs(c(b));
                    results.(sprintf('accuracy_%s_bin%d', flds{f}, b))(icnt)    = a(b);
                    
                    % repetition criterion as a function of pupil stuff
                    results.(sprintf('repetitioncrit_%s_bin%d', flds{f}, b))(icnt)  = -r(b);
                    
                    % repetition as a function of pupil stuff
                    results.(sprintf('repetition_%s_bin%d', flds{f}, b))(icnt)      = rep(b);
                    results.(sprintf('stimrepetition_%s_bin%d', flds{f}, b))(icnt)      = stimrep(b);
                    
                    % criterionshift
                    results.(sprintf('criterionshift_%s_bin%d', flds{f}, b))(icnt)  = cs(b);
                    
                    % next trial dprime and criterion
                    results.(sprintf('nextdprime_%s_bin%d', flds{f}, b))(icnt)      = nextd(b);
                    results.(sprintf('nextabscriterion_%s_bin%d', flds{f}, b))(icnt)   = abs(nextc(b));
                    
                end
            end
        end
    end
end

% correct repition for stimulus repetition by subtraction as well
results.repetitioncorr = results.repetition - results.stimrepetition;

end

function out = corrFunc(in)
% correlate two things, if there are too few datapoints output NaN

out = fisherz(corr(transpose(1:length(in)), in, ...
    'type', 'spearman', 'rows', 'complete'));
if abs(out) > 2, out = NaN; end

end
