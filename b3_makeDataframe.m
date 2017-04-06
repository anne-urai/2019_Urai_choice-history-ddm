% get a huge list with values for each participant
% can then work with this dataframe

clear; close all; clc;
addpath(genpath('~/code/Tools'));
warning off;

% ============================================ %
% TWO DIFFERENT DATASETS
% ============================================ %

usr = getenv('USER');
switch usr
    case 'anne' % local
        datasets = {'RT_RDK', 'projects/0/neurodec/Data/MEG-PL'};
    case 'aeurai' % lisa/cartesius
        datasets = {'RT_RDK', 'MEG-PL'};
end

set(groot, 'defaultaxesfontsize', 7, 'defaultaxestitlefontsizemultiplier', 1, ...
    'defaultaxestitlefontweight', 'bold', ...
    'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters');

for d = 1:length(datasets),
    
    % load data
    csvfile = dir(sprintf('~/Data/%s/HDDM/*allsj.csv', datasets{d}));
    alldata = readtable(sprintf('~/Data/%s/HDDM/%s', datasets{d}, csvfile.name));
    
    % rename some things
    try
        alldata.Properties.VariableNames{'trialnr'}     = 'trial';
        alldata.Properties.VariableNames{'blocknr'}     = 'block';
        
        alldata.session = alldata.session + 1; % start at 1
        % session 0 means the average of all sessions
    catch
        alldata.Properties.VariableNames{'subjnr'}   = 'subj_idx';
        alldata.Properties.VariableNames{'stim'}     = 'stimulus';
        alldata.Properties.VariableNames{'resp'}     = 'response';
        alldata.prevrt = circshift(alldata.rt, 1);
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
    
    % manually recode the drift rate parameters to match the specific session
    switch d
        case 1 % RT-RDK dataset
            sessions = 0:4;
        case 2
            sessions = 0:1;
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
            tab.(newvar)(tab.session == s + 1) = tab.(vars{v})(tab.session == 0);
            % remove the old one
            tab(:,{vars{v}}) = [];
        end
    end
    
    writetable(tab, sprintf('~/Data/%s/HDDM/summary/allindividualresults.csv', datasets{d}));
end