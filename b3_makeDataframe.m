% get a huge list with values for each participant
% can then work with this dataframe
function b3_makeDataframe(ds)

  close all; clc;
  addpath(genpath('~/code/Tools'));
  warning off;

  % ============================================ %
  % SUMMARIZE EACH DATASET
  % ============================================ %

  usr = getenv('USER');
  switch usr
  case 'anne' % local
    mypath = '~/Data/HDDM';
  case 'aeurai' % lisa/cartesius
    mypath = '/nfs/aeurai/HDDM';
  end
  
 datasets = {'Murphy', 'JW_yesno', 'NatComm', 'MEG', 'JW_PNAS', 'JW_fMRI', ...
  'Anke_2afc_sequential', 'Anke_MEG', 'Bharath_fMRI', 'Anke_merged'};
  datasets = {'MEG_MEGdata'};
 
  set(groot, 'defaultaxesfontsize', 7, 'defaultaxestitlefontsizemultiplier', 1, ...
  'defaultaxestitlefontweight', 'bold', ...
  'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters');

if ~exist('ds', 'var'), ds = 1:length(datasets); end

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
    results     = b3b_behaviouralMetrics(alldata);

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
    switch datasets{d}
    case 'RT_RDK'
      sessions = 1:5;
    case {'MEG', 'MEG_MEGsessions'};
      sessions = 1:5;
    case {'Anke_2afc_serial', 'Anke_2afc_neutral', 'Anke_2afc_repetitive', 'Anke_2afc_altenating'},
      sessions = 1:6;
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
        % remove the old one
        % tab(:,{vars{v}}) = [];
      end
    end

    % remove sessions where no data was recorded
    skippedSession = (isnan(nanmean(tab{:, 3:11}, 2)));
    tab(skippedSession, :) = [];

	switch whichFit
	case 1
    	writetable(tab, sprintf('%s/summary/%s/allindividualresults.csv', mypath, datasets{d}));
	case 2
    	writetable(tab, sprintf('%s/summary/%s/allindividualresults_Gsq.csv', mypath, datasets{d}));
	end
	
    fprintf('%s/summary/%s/allindividualresults.csv \n', mypath,  datasets{d});

  end
end
