function results = A5a_behaviouralMetrics()
% error vs correct, dprime vs criterion, serial choice bias
% as a function of pupil and RT

if ~isdeployed,
    addpath(genpath('~/code/MEG'));
    addpath(genpath('~/code/Tools'));
    addpath('~/Documents/fieldtrip');
    ft_defaults;
    addpath /Users/anne/Dropbox/code/MEG;
end
% warning('error', 'MATLAB:table:ModifiedVarnames');

% ========================================== %
% COMPUTE FOR ALL SUBJECTS
% ========================================== %

clear all; clc; close all;
subjectdata = subjectspecifics('ga');

% preallocate variables
varnames = {'subjnr', 'session', 'dprime', 'accuracy', 'criterion', 'abscriterion', 'rt', ...
    'pupil_correct', 'pupil_error', 'rt_correct', 'rt_error'};

results = array2table(nan(length(subjectdata.all)*3, length(varnames)), 'variablenames', varnames);
results.drug = repmat({'NaN'}, length(subjectdata.all)*3, 1);

metrics = {'dprime', 'criterion', 'abscriterion', 'accuracy', 'repetition', ...
    'stimrepetition', 'repetitioncrit', 'criterionshift', 'handshift', 'nextdprime', 'nextabscriterion'};

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

% preallocate
subjects = subjectdata.all;

% get all data
alldata = readtable(sprintf('%s/Data/CSV/2ifc_megdata_allsj.csv', subjectdata.path));
alldata.subjnr = round(alldata.idx / 1000000);
alldata.startHand(alldata.startHand > 20) = nan;

% for criterion shift
alldata.nextstim = circshift(alldata.stim, -1);
alldata.nextresp = circshift(alldata.resp, -1);
alldata.nextstim((diff(alldata.trial) ~= 1)) = NaN;
alldata.nextresp((diff(alldata.trial) ~= 1)) = NaN;

icnt = 0;
for sj = subjects,
    sjdat = subjectspecifics(sj);
    for s = 0:2,
        
        icnt                    = icnt + 1;
        results.subjnr(icnt)    = sj;
        results.session(icnt)   = s;
        
        switch s
            case 0
                % both sessions together
                data                    = alldata(alldata.subjnr == sj, :);
            otherwise
                data                    = alldata(alldata.subjnr == sj & alldata.session == s, :);
        end
        data(isnan(data.resp), :) = [];
        data.pupil              = data.decision_pupil;
        
        % ========================================== %
        % GENERAL STUFF
        % ========================================== %
        
        results.drug(icnt)          = {sjdat.drug};
        
        [d, c] = dprime(data.stim, data.resp);
        results.dprime(icnt)        = d;
        results.criterion(icnt)     = c;
        results.abscriterion(icnt)  = abs(c);
        results.accuracy(icnt)      = nanmean(data.correct);
        results.rt(icnt)            = nanmedian(data.rt);
        
        % some people don't have pupil data in each session
        if isempty(data),
            fprintf('skipping sj %d, session %d \n', sj, s);
            continue;
        end
        
        % measure of repetition behaviour
        data.repeat = [~(abs(diff(data.resp)) > 0); NaN];
        data.stimrepeat = [~(abs(diff(data.stim)) > 0); NaN];
        
        % skip trials at boundaries
        data.repeat((diff(data.trial) ~= 1)) = NaN;
        data.stimrepeat((diff(data.trial) ~= 1)) = NaN;
        
        results.repetition(icnt)        = nanmean(data.repeat);
        results.stimrepetition(icnt)    = nanmean(data.stimrepeat);
        
        % criterion based on repetition and stimulus sequences
        [~, c] = dprime(data.stimrepeat, data.repeat);
        results.repetitioncrit(icnt)    = -c;
        
        % criterion based on next trial bias, then collapsed
        results.criterionshift(icnt)    = criterionshift(data.resp, data.nextstim, data.nextresp);
        % does the random hand they press cause a bias?
        results.handshift(icnt)         = criterionshift(data.startHand, data.stim, data.resp);
        
        % pupilstuff
        results.pupil_error(icnt)       = nanmean(data.pupil(data.correct == 0));
        results.pupil_correct(icnt)     = nanmean(data.pupil(data.correct == 1));
        
        results.rt_error(icnt)          = nanmedian(data.rt(data.correct == 0));
        results.rt_correct(icnt)        = nanmedian(data.rt(data.correct == 1));
        
        % ========================================== %
        % DPRIME, CRITERION, ABSOLUTE CRITERION
        % correlation to pupil and RT
        % ========================================== %
        
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
            [d, c]      = splitapply(@dprime, data.stim, data.resp, binIdx);
            a           = splitapply(@nanmean, data.correct, binIdx);
            [~, r]      = splitapply(@dprime, data.stimrepeat, data.repeat, binIdx);
            rep         = splitapply(@nanmean, data.repeat, binIdx);
            stimrep     = splitapply(@nanmean, data.stimrepeat, binIdx);
            cs          = splitapply(@criterionshift, data.resp, data.nextstim, data.nextresp, binIdx);
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
            results.([flds{f} '_nextabscriterion_corr'])(icnt)     = corrFunc(abs(nextc));
            
            % ========================================== %
            % divide into 3 bins for viz
            % ========================================== %
            
            binIdx = discretize(useFld, [-inf quantile(useFld, 2) inf]);
            
            [d, c] = splitapply(@dprime, data.stim, data.resp, binIdx);
            a      = splitapply(@nanmean, data.correct, binIdx);
            [~, r] = splitapply(@dprime, data.stimrepeat, data.repeat, binIdx);
            rep    = splitapply(@nanmean, data.repeat, binIdx);
            cs     = splitapply(@criterionshift, data.resp, data.nextstim, data.nextresp, binIdx);
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

% correct repition for stimulus repetition by subtraction as well
results.repetitioncorr = results.repetition - results.stimrepetition;

end

function out = corrFunc(in)
% correlate two things, if there are too few datapoints output NaN

out = fisherz(corr(transpose(1:length(in)), in, ...
    'type', 'spearman', 'rows', 'complete'));
if abs(out) > 2, out = NaN; end

end


