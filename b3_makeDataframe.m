% get a huge list with values for each participant
% can then work with this dataframe

clear; close all; clc;
addpath(genpath('~/code/Tools'));
warning off;

% ============================================ %
% SUMMARIZE EACH DATASET
% ============================================ %

usr = getenv('USER');
switch usr
    case 'anne' % local
        datasets = {'RT_RDK/HDDM', 'projects/0/neurodec/Data/MEG-PL/HDDM', ...
            'projects/0/neurodec/Data/MEG-PL/HDDM-S1', 'projects/0/neurodec/Data/MEG-PL/HDDM-S2', ...
            'Anke_2afc_sequential/HDDM-alternating', 'Anke_2afc_sequential/HDDM-neutral', ...
            'Anke_2afc_sequential/HDDM-repetitive'};
    case 'aeurai' % lisa/cartesius
        datasets = {'RT_RDK', 'MEG-PL'};
end

set(groot, 'defaultaxesfontsize', 7, 'defaultaxestitlefontsizemultiplier', 1, ...
    'defaultaxestitlefontweight', 'bold', ...
    'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters');

for d = 2; %3:length(datasets),
    
    % load data
    csvfile = dir(sprintf('~/Data/%s/*.csv', datasets{d}));
    alldata = readtable(sprintf('~/Data/%s/%s', datasets{d}, csvfile.name));
    
    % rename some things
    try
        alldata.Properties.VariableNames{'trialnr'}     = 'trial';
        alldata.Properties.VariableNames{'blocknr'}     = 'block';
        
        alldata.session = alldata.session + 1; % start at 1
        % session 0 means the average of all sessions
    catch
        try
            alldata.Properties.VariableNames{'subjnr'}   = 'subj_idx';
            alldata.Properties.VariableNames{'stim'}     = 'stimulus';
            alldata.Properties.VariableNames{'resp'}     = 'response';
            alldata.prevrt = circshift(alldata.rt, 1);
        end
    end
    
    % compute a bunch of basic things from Matlab
    results     = b3b_behaviouralMetrics(alldata);
    
    % get the summary results from HDDM
    hddmresults = readtable(sprintf('~/Data/%s/summary/individualresults.csv', datasets{d}));
    
    % most parameters will go under session 0
    hddmresults.session = zeros(size(hddmresults.subjnr));
    
    % will only keep session 0 stuff
    allresults = innerjoin(results, hddmresults);
    
    % now add back all the stuff from the different sessions
    allresults2 = tableAppend(allresults, results);
    
    % remove duplicate rows, save only those with HDDM info
    % http://stackoverflow.com/questions/27547463/matlab-delete-duplicate-table-entries-on-multiple-columns
    [~, ind] = unique(allresults2(:, [1 2]), 'rows');
    tab      = allresults2(ind,:);
    
    % manually recode the drift rate parameters to match the specific session
    switch d
        case 1 % RT-RDK dataset
            sessions = 0:4;
        case 2
            sessions = 0:1;
        otherwise
            % for Anke's data, skip this
            sessions = [];
    end
    
    for s = sessions,
        vars = find(~cellfun(@isempty, strfind(tab.Properties.VariableNames, sprintf('v_%d', s))));
        vars = tab.Properties.VariableNames(vars);
        
        for v = 1:length(vars),
            
            % if this is the first session, make a new column for
            % the overall drift rate (which will then be repopulated per
            % session)
            if s == min(sessions),
                newvar = regexprep(vars{v}, sprintf('_%d__', s), '__');
                tab.(newvar) = nan(size(tab.(vars{v})));
            end
            
            % then, move the values over
            tab.(newvar)(tab.session == s) = tab.(vars{v})(tab.session == 0);
            % remove the old one
            % tab(:,{vars{v}}) = [];
        end
    end
    
    % remove sessions where no data was recorded
    skippedSession = (isnan(nanmean(tab{:, 3:11}, 2)));
    tab(skippedSession, :) = [];
    
    writetable(tab, sprintf('~/Data/%s/summary/allindividualresults.csv', datasets{d}));
    
end

% ============================================ %
% APPEND DATASETS THAT GO TOGETHER
% ============================================ %

clearvars -except datasets

% 1. MEG-PL HDDM, just stimcoding
dat1 = readtable(sprintf('~/Data/%s/summary/allindividualresults.csv', datasets{3}));
dat1.session = dat1.session + 1;
dat2 = readtable(sprintf('~/Data/%s/summary/allindividualresults.csv', datasets{4}));
dat2(dat2.session == 1, :) = [];
dat2.session = dat2.session + 2;
dat2.Properties.VariableNames{'v_1__stimcoding_prevresp_dc_z'} = 'v_0__stimcoding_prevresp_dc_z';
dat2 = dat2(:, dat1.Properties.VariableNames);
dat3 = cat(1, dat1, dat2);
writetable(dat3, sprintf('~/Data/%s/summary/allindividualresults_separatesessions.csv', datasets{2}));

% 2. Anke's data
clearvars -except datasets
dat1 = readtable(sprintf('~/Data/%s/summary/allindividualresults.csv', datasets{5}));
dat1.transprob = 0.2 * ones(size(dat1.subjnr)); % alternating
dat2 = readtable(sprintf('~/Data/%s/summary/allindividualresults.csv', datasets{6}));
dat2.transprob = 0.5 * ones(size(dat2.subjnr)); % neutral
dat3 = readtable(sprintf('~/Data/%s/summary/allindividualresults.csv', datasets{7}));
dat3.transprob = 0.8 * ones(size(dat3.subjnr)); % repetitive
dat = cat(1, dat1, dat2, dat3);
writetable(dat, sprintf('~/Data/%s/allindividualresults.csv', 'Anke_2afc_sequential/HDDM'));
