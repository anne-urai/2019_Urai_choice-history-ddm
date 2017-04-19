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
        datasets = {'RT_RDK', 'projects/0/neurodec/Data/MEG-PL', 'Anke_2afc_sequential'};
    case 'aeurai' % lisa/cartesius
        datasets = {'RT_RDK', 'MEG', 'Anke_serial'};
end

set(groot, 'defaultaxesfontsize', 7, 'defaultaxestitlefontsizemultiplier', 1, ...
    'defaultaxestitlefontweight', 'bold', ...
    'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters');

for d = 3:length(datasets),
    disp(datasets{d});
    
    % load data
    csvfile = dir(sprintf('~/Data/%s/HDDM/*.csv', datasets{d}));
    csvfile = csvfile(arrayfun(@(x) ~strcmp(x.name(1),'.'), csvfile)); % remove hidden files
    alldata = readtable(sprintf('~/Data/%s/HDDM/%s', datasets{d}, csvfile.name));
    
    % recode Anke's stimulus into stim and coh
    if d == 3,
        alldata.coherence   = abs(alldata.stimulus);
        alldata.stimulus2   = sign(alldata.stimulus);
        alldata.stimulus2(alldata.coherence == 0) = sign(alldata.motionenergy(alldata.coherence == 0));
        alldata.stimulus    = alldata.stimulus2;
    end
    
    % compute a bunch of basic things from Matlab
    results     = b3b_behaviouralMetrics(alldata);
    
    % get the summary results from HDDM
    hddmresults = readtable(sprintf('~/Data/%s/HDDM/summary/individualresults.csv', datasets{d}));
    
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
    
    % ============================================ %
    % RECODE SESSION-SPECIFIC PARAMETERS
    % ============================================ %
    
    % manually recode the drift rate parameters to match the specific session
    switch d
        case 1 % RT-RDK dataset
            sessions = 1:5;
        case 2 % MEG
            sessions = 1:5;
        case 3 % Anke
            sessions = 1:6;
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
            % remove the old one
            % tab(:,{vars{v}}) = [];
        end
    end
    
    % remove sessions where no data was recorded
    skippedSession = (isnan(nanmean(tab{:, 3:11}, 2)));
    tab(skippedSession, :) = [];
    
    writetable(tab, sprintf('~/Data/%s/HDDM/summary/allindividualresults.csv', datasets{d}));
    
end
