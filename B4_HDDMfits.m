% get the data from two HDDM models that were run

clear all; close all;
usepath = '/Users/anne/Data/projects/0/neurodec/Data/MEG-PL/HDDM';
mdls = {'stimcoding', 'stimcoding_prevresp_z', 'stimcoding_prevresp_dc', ...
    'regress_dc', 'regress_dc_prevresp', 'regress_dc_prevresp_prevpupil_prevrt'};

for m = 5:length(mdls),
    
    % ============================================ %
    % CHECK MODEL CONVERGENCE
    % ============================================ %
    
    % 1. Gelman-Rubin rhat
    % Values should be close to 1 and not larger than 1.02 which would indicate convergence problems.
    %     gr = readtable(sprintf('%s/%s/gelman_rubin.txt', usepath, mdls{m}), ...
    %         'delimiter', ':', 'readvariablenames', false);
    %     if ~all(gr{:, 2} < 1.02) || ~all(gr{:, 2} > 0.98)
    %         warning('rhat indicates non-convergence');
    %     end
    %
    % ============================================ %
    % GET PARAMETERS OF INTEREST
    % ============================================ %
    
    % alles tussen haakjes () zijn condities (session/prevresp)
    % alles
    
    res  = readtable(sprintf('%s/%s/results-combined.csv', usepath, mdls{m}));
    varn = res{:, 1};
    
    % indicate the parameter for each row
    parameter      = regexp(varn, '^[a-z]+(?=[\.\(\_])', 'match');
    parameter(cellfun(@isempty, parameter)) = {'NaN'};
    parameter      = [parameter{:}]';
    res.parameter  = parameter;
    
    % parse the rownames to construct new variables
    subjnr      = regexp(varn, '(?<=\.)[0-9.]*', 'match');
    subjnr(cellfun(@isempty, subjnr)) = {'NaN'};
    for s = 1:length(subjnr),
        if iscell(subjnr{s}) && numel(subjnr{s}) > 1,
            if strcmp(subjnr{s}{2}, '0.0'),
                subjnr{s} = subjnr{s}{1};
            end
        end
    end
    
    subjnr      = [subjnr{:}]';
    res.subjnr  = cellfun(@str2num, subjnr, 'un', 1);
    
    % also add the session nrs for the variable that was additionally split by previous response
    sessionnr    = regexp(varn, '(?<=\().*(?=\))', 'match');
    sessionnr(cellfun(@isempty, sessionnr)) = {'NaN'};
    sessionnr    = [sessionnr{:}]';
    res.session  = cellfun(@str2num, sessionnr, 'un', 1);
    
    %     sessionnr    = regexp(varn, '(?<=0\.)[0-9-]*(?=\.0)', 'match');
    %     idx          = find(~cellfun(@isempty, sessionnr));
    %     sessionnr(cellfun(@isempty, sessionnr)) = {'NaN'};
    %     sessionnr    = [sessionnr{:}]';
    %     res.session(idx) = cellfun(@str2num, sessionnr(idx), 'un', 1);
    %
    %     % code for the previous response
    %     prevresp    = regexp(varn, '(?<=\()[0-9-]*(?=\.)', 'match');
    %     prevresp(cellfun(@isempty, prevresp)) = {'0'};
    %     prevresp    = [prevresp{:}]';
    %     res.prevresp  = cellfun(@str2num, prevresp, 'un', 1);
    
    % remove some crap
    res(:, {'std', 'x2_5q', 'x25q', 'x50q', 'x75q', 'x97_5q', 'mcErr', 'Var1'}) = [];
    
    % change into a wide format
    res2 = unstack(res, 'mean', 'parameter');
    res2(:, {'NaN'}) = [];
    writetable(res2, sprintf('%s/%s/results-parsed.csv', usepath, mdls{m}));
    
end

% ============================================ %
% PLOT DEPENDENCE ON PREVIOUS TRIAL RESPONSE
%% ============================================ %

clearvars -except usepath mdls
cnt = 1;
for m = 1:length(mdls),
    res = readtable(sprintf('%s/%s/results-parsed.csv', usepath, mdls{m}));
    
    switch m
        case 1
            par = 'dc';
        case 2
            par = 'z';
    end
    
    for s = 1:2,
        prevresp1 = res.(par)(res.session == s & res.prevresp == 1);
        prevresp0 = res.(par)(res.session == s & res.prevresp == -1);
        
        subplot(4,4,cnt); cnt = cnt + 1;
        plot([prevresp0 prevresp1]');
        pval = permtest(prevresp0, prevresp1);
        mysigstar(gca, [1 2], max(get(gca, 'ylim')), pval);
        ylabel(par);
        xlabel('Previous response');
        set(gca, 'xtick', 1:2, 'xticklabel', [-1 1], 'xlim', [0.5 2.5]);
        box off; title(sprintf('Session %d', s));
    end
end

% correlate the previous response's effect on z and dc

% correlate both with pure repetition probability

