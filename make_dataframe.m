
function make_dataframe(datasets)

% Code to fit the history-dependent drift diffusion models described in
% Urai AE, Gee JW de, Donner TH (2018) Choice history biases subsequent evidence accumulation. bioRxiv:251595
%
% MIT License
% Copyright (c) Anne Urai, 2018
% anne.urai@gmail.com

close all; clc;
addpath(genpath('~/code/Tools'));
warning off; global mypath

% ============================================ %
% SUMMARIZE EACH DATASET
% ============================================ %

set(groot, 'defaultaxesfontsize', 7, 'defaultaxestitlefontsizemultiplier', 1, ...
    'defaultaxestitlefontweight', 'bold', ...
    'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters');

ds = 1:length(datasets);

for d = ds,
    disp(datasets{d});

    % load data
    csvfile = dir(sprintf('%s/%s/*.csv', mypath, datasets{d}));
    csvfile = csvfile(arrayfun(@(x) ~strcmp(x.name(1),'.'), csvfile)); % remove hidden files
    alldata = readtable(sprintf('%s/%s/%s', mypath, datasets{d}, csvfile.name));

    % recode Anke's stimulus into stim and coh
    if ~isempty(strfind(datasets{d}, 'Anke')) | ~isempty(strfind(datasets{d}, 'NatComm')) | ~isempty(strfind(datasets{d}, 'Bharath')),
        alldata.coherence   = abs(alldata.stimulus);
        alldata.stimulus2   = sign(alldata.stimulus);
        try; alldata.stimulus2(alldata.coherence == 0) = sign(alldata.motionenergy(alldata.coherence == 0)); end
        alldata.stimulus    = alldata.stimulus2;
    end

    % compute a bunch of basic things from Matlab
    results     = define_behavioral_metrics(alldata);
    results     = [results; array2table(nan(size(results.Properties.VariableNames)), ...
        'variablenames', results.Properties.VariableNames)];
    results.session(isnan(results.subjnr)) = 0;
    results.repetition(isnan(results.subjnr)) = 0;
    results.subjnr(isnan(results.subjnr)) = 0;

    % add personality scores and drug conditions
    switch datasets{d}
        case {'MEG', 'MEG_MEGsessions', 'MEG_MEGdata'},
            disp('adding in personality questionnaires');
            results.drug = repmat({'NaN'}, length(results.dprime), 1);
            results.BIS = nan(size(results.dprime));
            results.BAS = nan(size(results.dprime));
            results.AQ = nan(size(results.dprime));
            results.PSWQ = nan(size(results.dprime));

            sjs = unique(results.subjnr)';
            sjs(sjs == 0) = []; % exclude group average
            for sj = sjs,
                subjectdata = subjectspecifics(sj);
                results.drug(results.subjnr == sj)  = {subjectdata.drug};
                results.BIS(results.subjnr == sj)   = subjectdata.BIS;
                results.BAS(results.subjnr == sj)   = subjectdata.BAS;
                results.AQ(results.subjnr == sj)    = subjectdata.AQ;
                results.PSWQ(results.subjnr == sj)  = subjectdata.PSWQ;
            end
    end

    for whichFit = 1:2,

        switch whichFit
            case 1
                % get the summary results from HDDM
                hddmresults = readtable(sprintf('%s/summary/%s/individualresults.csv', mypath, datasets{d}));
            case 2
                if exist(sprintf('%s/summary/%s/individualresults_Gsq.csv', mypath, datasets{d}), 'file'),
                    hddmresults = readtable(sprintf('%s/summary/%s/individualresults_Gsq.csv', mypath, datasets{d}));
                else
                    continue;
                end
        end

        % most parameters will go under session 0
        hddmresults.session = zeros(size(hddmresults.subjnr));

        % will only keep session 0 stuff
        allresults = innerjoin(results, hddmresults, 'keys', {'subjnr', 'session'});

        % now add back all the stuff from the different sessions
        allresults2 = tableAppend(allresults, results);

        % remove duplicate rows, save only those with HDDM info
        % http://stackoverflow.com/questions/27547463/matlab-delete-duplicate-table-entries-on-multiple-columns
        [~, ind] = unique(allresults2(:, [1 2]), 'rows');
        tab      = allresults2(ind,:);
        assert(any(ismember(tab.subjnr, 0)), 'group average missing');

        % ============================================ %
        % RECODE SESSION-SPECIFIC PARAMETERS
        % ============================================ %

        % manually recode the drift rate parameters to match the specific session
        switch datasets{d}
            case 'RT_RDK'
                sessions = 1:5;
            case {'MEG', 'MEG_MEGsessions'};
            case 'NatComm'
                sessions = 1:5;
        end

        varidx = find(~cellfun(@isempty, strfind(tab.Properties.VariableNames, sprintf('_s%d_', 1))));
        vars   = tab.Properties.VariableNames(varidx);

        for v = 1:length(vars),
            for s = sessions,

                % if this is the first session, make a new column for
                % the overall drift rate (which will then be repopulated per
                % session)
                if s == min(sessions),
                    newvar = regexprep(vars{v}, sprintf('_s%d__', s), '__');
                    tab.(newvar) = nan(size(tab.(vars{v})));
                    thisvar = vars{v};
                else
                    thisvar = regexprep(vars{v}, '_s1_', sprintf('_s%d_', s));
                end

                % then, move the values over
                try
                    tab.(newvar)(tab.session == s) = tab.(thisvar)(tab.session == 0);
                    % can happen that there is no session 2-4 (MEG pupil)
                end
            end
        end

        % remove sessions where no data was recorded
        skippedSession = (isnan(nanmean(tab{:, 3:11}, 2)));
        tab(skippedSession, :) = [];

        % group-level HDDM estimates - remove sjnr and session nr
        tab.session(tab.subjnr == 0) = NaN;
        tab.subjnr(tab.subjnr == 0) = NaN;

        % ============================================ %
        % SAVE TO FIGSHARE FOR CLEARER OVERVIEW
        % ============================================ %

        switch whichFit
            case 1
                writetable(tab, sprintf('%s/summary/%s/allindividualresults.csv', mypath, datasets{d}));

                %% ALSO SAVE FOR FIGSHARE
                tab2 = tab(tab.session == 0 | isnan(tab.session), :);

                switch d
                case 1
                    writetable(tab2, sprintf('%s/summary/visual_motion_2afc_rt_hddmfits.csv', mypath));
                case 2
                    writetable(tab2, sprintf('%s/summary/visual_motion_2afc_fd_hddmfits.csv', mypath));
                case 3
                    writetable(tab2, sprintf('%s/summary/visual_motion_2ifc_fd_1_hddmfits.csv', mypath));
                case 4
                    writetable(tab2, sprintf('%s/summary/visual_motion_2ifc_fd_2_hddmfits.csv', mypath));
                case 5
                    writetable(tab2, sprintf('%s/summary/visual_contrast_yesno_hddmfits.csv', mypath));
                case 6
                    writetable(tab2, sprintf('%s/summary/auditory_yesno_hddmfits.csv', mypath));
                end

            case 2
                writetable(tab, sprintf('%s/summary/%s/allindividualresults_Gsq.csv', mypath, datasets{d}));

                %% ALSO SAVE FOR FIGSHARE
                tab2 = tab(tab.session == 0, :);

                switch d
                case 1
                    writetable(tab2, sprintf('%s/summary/visual_motion_2afc_rt_gsquarefits.csv', mypath));
                case 2
                    writetable(tab2, sprintf('%s/summary/visual_motion_2afc_fd_gsquarefits.csv', mypath));
                case 3
                    writetable(tab2, sprintf('%s/summary/visual_motion_2ifc_fd_1_gsquarefits.csv', mypath));
                case 4
                    writetable(tab2, sprintf('%s/summary/visual_motion_2ifc_fd_2_gsquarefits.csv', mypath));
                case 5
                    writetable(tab2, sprintf('%s/summary/visual_contrast_yesno_gsquarefits.csv', mypath));
                case 6
                    writetable(tab2, sprintf('%s/summary/auditory_yesno_gsquarefits.csv', mypath));
                end
        end
        fprintf('%s/summary/%s/allindividualresults.csv \n', mypath,  datasets{d});



    end
end
