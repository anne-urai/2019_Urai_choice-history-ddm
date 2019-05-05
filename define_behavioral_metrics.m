function results = define_behavioral_metrics(alldata)

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
varnames = {'subjnr', 'session', 'dprime', 'accuracy', 'criterion', 'repetition', 'repetition2', 'repetition3', ...
    'criterionshift', 'repetition_prevcorrect', 'repetition_preverror', 'posterrorslowing', ...
    'repetition_congruent', 'repetition_incongruent', 'repetition_fastRT', 'repetition_slowRT'};

nrSess          = length(unique(alldata.session)) + 1;
results         = array2table(nan(length(unique(alldata.subj_idx))*nrSess, length(varnames)), 'variablenames', varnames);

% preallocate dprime for different coherence levels
if sum(strcmp(alldata.Properties.VariableNames, 'coherence')) > 0,
    cohlevels = unique(alldata.coherence);
    for c = 1:length(cohlevels),
        vrnm = ['dprime_c' num2str(cohlevels(c)*100)];
        vrnm = regexprep(vrnm, '\.', '\_'); % replace points in varname
        results.(vrnm) = nan(size(results.dprime));
        
        % ALSO COMPUTE REPETITION SEPARATELY FOR EACH LEVEL OF COHERENCE
        vrnm = ['repetition_c' num2str(cohlevels(c)*100)];
        vrnm = regexprep(vrnm, '\.', '\_'); % replace points in varname
        results.(vrnm) = nan(size(results.dprime));
        
    end
end

% % get all data
subjects      = unique(alldata.subj_idx)';

% recode correct
if all(cellfun(@isempty, strfind(alldata.Properties.VariableNames, 'correct'))),
    tmpstim = (alldata.stimulus > 0);
    tmpresp = (alldata.response > 0);
    alldata.correct = (tmpstim == tmpresp);
end

% for mulder et al. analysis
alldata.prevstim = circshift(alldata.stimulus, 1);
alldata.prevresp = circshift(alldata.response, 1);
alldata.congruent = (sign(alldata.prevresp - 0.1) == alldata.stimulus);

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
   % for s = [0 unique(alldata.session)'],
    for s = 0,
        
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
        % results.abscriterion(icnt)  = abs(c);
        results.accuracy(icnt)      = nanmean(data.correct);
        results.rt(icnt)            = nanmedian(data.rt);
        results.bias(icnt)          = nanmean(data.response);
        
        % post-error slowing
        prevcorrect = (data.prevstim == data.prevresp);
        results.posterrorslowing(icnt) = nanmean(data.rt(prevcorrect == 0)) - ...
            nanmean(data.rt(prevcorrect == 1));
        
        % measure of repetition behaviour
        % data.repeat = [~(abs(diff(data.response)) > 0); NaN];
        % data.stimrepeat = [~(abs(diff(data.stimulus)) > 0); NaN];
        
        % 01.10.2017, use the same metric as in MEG, A1c_writeCSV.m
        for l = 1:6,
            data.(['prev' num2str(l) 'resp']) = circshift(data.response, l);
            data.(['prev' num2str(l) 'stim']) = circshift(data.stimulus, l);
            data.(['prev' num2str(l) 'resp'])(data.(['prev' num2str(l) 'resp']) == 0) = -1;

            data.(['repeat' num2str(l)]) = double(data.response == circshift(data.response, l));
            
            % EXCLUDE TRIALS THAT ARE NOT CONTINUOUS
            wrongTrls = ((data.trial - circshift(data.trial, l)) ~= l);
            data.(['repeat' num2str(l)])(wrongTrls) = NaN;
        end
        
        data.repeat = data.repeat1;
        data.stimrepeat = [NaN; (diff(data.response) == 0)];
       
        if sum(strcmp(data.Properties.VariableNames, 'coherence')) > 0,
            cohlevels = unique(data.coherence);
            for c = 1:length(cohlevels),
                vrnm = ['dprime_c' num2str(cohlevels(c)*100)];
                vrnm = regexprep(vrnm, '\.', '\_'); % replace points in varname
                % disp(vrnm);
                results.(vrnm)(icnt) = ...
                    dprime(data.stimulus(data.coherence == cohlevels(c)), ...
                    data.response(data.coherence == cohlevels(c)));
                if cohlevels(c) > 0,
                    assert(~isnan(results.(vrnm)(icnt)));
                end
                
                %% ALSO COMPUTE REPETITION SEPARATELY FOR EACH LEVEL OF COHERENCE
                vrnm = ['repetition_c' num2str(cohlevels(c)*100)];
                vrnm = regexprep(vrnm, '\.', '\_'); % replace points in varname
                results.(vrnm)(icnt) = nanmean(data.repeat(data.coherence == cohlevels(c)));
            end
        end
        
        % ======================================= %
        % add repetition across longer lags
        % for figure 6c
        % ======================================= %
        
        for l = 1:6,
            results.(['repetition' num2str(l)])(icnt) = nanmean(data.(['repeat' num2str(l)]));
            results.(['repetition_correct' num2str(l)])(icnt) = ...
            nanmean(data.(['repeat' num2str(l)])((data.(['prev' num2str(l) 'stim']) > 0) == (data.(['prev' num2str(l) 'resp']) > 0)));
            results.(['repetition_error' num2str(l)])(icnt) = ...
            nanmean(data.(['repeat' num2str(l)])((data.(['prev' num2str(l) 'stim']) > 0) ~= (data.(['prev' num2str(l) 'resp']) > 0)));
        end 

        % for the first lag, no correction (perhaps not necessary?)
        results.repetition(icnt)        = nanmean(data.repeat1);
        
        % also compute this after error and correct trials
        results.repetition_prevcorrect(icnt) = nanmean(data.repeat((data.prevstim > 0) == (data.prevresp > 0)));
        results.repetition_preverror(icnt)   = nanmean(data.repeat((data.prevstim > 0) ~= (data.prevresp > 0)));

        % split by congruency
        results.repetition_congruent(icnt) = nanmean(data.repeat(data.congruent == 1));
        results.repetition_incongruent(icnt) = nanmean(data.repeat(data.congruent == 0));
        
        % REPETITION FOR SLOW VS FAST RTS
        results.repetition_fastRT(icnt) = nanmean(data.repeat(data.rt < nanmedian(data.rt)));
        results.repetition_slowRT(icnt) = nanmean(data.repeat(data.rt > nanmedian(data.rt)));

        % get the cumulative p(repeat) after a sequence, as in hermoso-mendizabal fig 2e
        results.cum_rep0(icnt) = nanmean(data.repeat);
        for sl = 1:6,
            startIndex = strfind(data.repeat', ones(1,sl));
            endIndex = startIndex + sl;
            endIndex(endIndex > length(data.repeat)) = [];
            try
                assert(numel(endIndex) > 5, 'not enough sequences to average over');
                results.(['cum_rep' num2str(sl)])(icnt) = nanmean(data.repeat(endIndex));
                results.(['cum_rep' num2str(sl)])(icnt) = NaN;
            end
        end
     
        % and alternating sequences
        results.cum_alt0(icnt) = nanmean(data.repeat);
        for sl = 1:6,
            startIndex = strfind(data.repeat', zeros(1,sl));
            endIndex = startIndex + sl;
            endIndex(endIndex > length(data.repeat)) = [];
            assert(numel(endIndex) > 5, 'not enough sequences to average over');
            results.(['cum_alt' num2str(sl)])(icnt) = nanmean(data.repeat(endIndex));
        end
    end
end

end

function out = corrFunc(in)
% correlate two things, if there are too few datapoints output NaN

out = fisherz(corr(transpose(1:length(in)), in, ...
    'type', 'spearman', 'rows', 'complete'));
if abs(out) > 2, out = NaN; end

end
