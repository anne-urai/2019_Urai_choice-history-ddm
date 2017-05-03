function e6_serialBias_SfN_modelFree_CRF
% ========================================== %
% conditional response functions from White & Poldrack
% ========================================== %

addpath(genpath('~/code/Tools'));
warning off; close all; clear;
global datasets datasetnames
datasets = {'RT_RDK', 'projects/0/neurodec/Data/MEG-PL', 'Anke_2afc_neutral', 'NatComm'};
datasetnames = {'RT', '2IFC', 'Anke neutral', 'NatComm'};

set(groot, 'defaultaxesfontsize', 6, 'defaultaxestitlefontsizemultiplier', 1.1, ...
    'defaultaxestitlefontweight', 'bold', ...
    'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters');

% ========================================== %
% MODELFREE MEASURE OF BIAS
% RT DISTRIBUTIONS
% ========================================== %

clf;
for d = 1:length(datasets),
    
    % load data
    clearvars -except d datasets cnt datasetnames whichSJ ylims;
    csvfile = dir(sprintf('~/Data/%s/HDDM/*.csv', datasets{d}));
    csvfile = csvfile(arrayfun(@(x) ~strcmp(x.name(1),'.'), csvfile)); % remove hidden files
    alldata = readtable(sprintf('~/Data/%s/HDDM/%s', datasets{d}, csvfile.name));
    
    % define repeaters and alternators based on dc
    dat = readtable(sprintf('~/Data/%s/HDDM/summary/allindividualresults.csv', ...
        datasets{d}));
    if d > 2,
        alldata.stimulus = sign(alldata.motionenergy);
    end
    
    % find alternators and repeaters
    allsubjects = unique(dat.subjnr); clear groups;
    groups{1}   = unique(dat.subjnr((dat.v_prevresp__regressdcprevrespstim < 0 & ...
        dat.criterionshift < 0)));
    groups{3}   = unique(dat.subjnr((dat.v_prevresp__regressdcprevrespstim > 0)  & ...
        dat.criterionshift > 0));
    groups{2}   = setdiff(allsubjects, [groups{1}; groups{3}]);
    alldata.session(:) = 1;
    
    cnt = d;
    groupnames  = {'alternators', 'rest', 'repeaters'};
    for gr = 1:length(groups),
        sph{gr} = subplot(4,4,cnt); cnt = cnt + 4;
        h = plotCRF(alldata(ismember(alldata.subj_idx, groups{gr}), :));
        title({datasetnames{d}, sprintf('%s n = %d', groupnames{gr}, numel(groups{gr}))});
        axis square;
    end
    
end

%suplabel(['Data: ' datasetnames{d}], 't', [sph{:}]);
l               = legend(h, {'Prev choice -1', 'Prev choice 1'}, 'Location', 'South');
l.Box           = 'off';
l.Position(2)   = l.Position(2) - 0.15;

print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/fig4_CRF.pdf'));


% ========================================== %
% MODELFREE MEASURE OF BIAS
% RT DISTRIBUTIONS
% ========================================== %
if 0,
    
    close;
    for d = 1:length(datasets),
        
        % load data
        clearvars -except d datasets cnt datasetnames whichSJ ylims;
        csvfile = dir(sprintf('~/Data/%s/HDDM/*.csv', datasets{d}));
        csvfile = csvfile(arrayfun(@(x) ~strcmp(x.name(1),'.'), csvfile)); % remove hidden files
        alldata = readtable(sprintf('~/Data/%s/HDDM/%s', datasets{d}, csvfile.name));
        
        % define repeaters and alternators based on dc
        dat = readtable(sprintf('~/Data/%s/HDDM/summary/allindividualresults.csv', ...
            datasets{d}));
        if d == 3,
            alldata.stimulus = sign(alldata.motionenergy);
        end
        % find alternators and repeaters
        groups{1}   = unique(dat.subjnr((dat.v_prevresp__regressdcprevrespstim < 0)));
        groups{2}   = unique(dat.subjnr((dat.v_prevresp__regressdcprevrespstim > 0)));
        
        cnt = d;
        groupnames  = {'alternators', 'repeaters'};
        for gr = 1:length(groups),
            sph{gr} = subplot(4,4,cnt); cnt = cnt + 4;
            h = plotRTdist(alldata(ismember(alldata.subj_idx, groups{gr}), :));
            title({datasetnames{d}, sprintf('%s n = %d', groupnames{gr}, numel(groups{gr}))});
            axis square;
        end
    end
    
    %suplabel(['Data: ' datasetnames{d}], 't', [sph{:}]);
    l               = legend(h, {'Prev choice -1', 'Prev choice 1'}, 'Location', 'South');
    l.Box           = 'off';
    l.Position(2)   = l.Position(2) - 0.15;
    
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/fig4_RTdist.pdf'));
    
    % ========================================== %
    % SINGLE SUBJECT PLOTS
    % ========================================== %
    
    for d = 1:length(datasets),
        
        % load data
        close all;
        clearvars -except d datasets cnt datasetnames whichSJ ylims;
        csvfile = dir(sprintf('~/Data/%s/HDDM/*.csv', datasets{d}));
        csvfile = csvfile(arrayfun(@(x) ~strcmp(x.name(1),'.'), csvfile)); % remove hidden files
        alldata = readtable(sprintf('~/Data/%s/HDDM/%s', datasets{d}, csvfile.name));
        
        % define repeaters and alternators based on dc
        dat = readtable(sprintf('~/Data/%s/HDDM/summary/allindividualresults.csv', ...
            datasets{d}));
        if d == 3,
            alldata.stimulus = sign(alldata.motionenergy);
        end
        
        % sort by serial choice bias
        nrsubpl     = ceil(sqrt(numel(unique(alldata.subj_idx))));
        [~, idx]    = sort(dat.v_prevresp__regressdcprevrespstim(dat.session == 0));
        sjs         = dat.subjnr(dat.session == 0);
        subjects    = sjs(idx);
        
        cnt = 1;
        for sj = subjects',
            subplot(nrsubpl, nrsubpl, cnt); cnt = cnt + 1;
            h = plotCRF(alldata(alldata.subj_idx == sj, :));
            title(sprintf('P%02d', sj));
            axis square;
        end
        print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/fig4_CRF_data%d_allsubjects.pdf', d));
    end
end

end

% ========================================== %
% function for conditional response function
% white & poldrack 2014
% ========================================== %

function h = plotCRF(data)

hold on;
qntls = [0 0.1 0.3 0.5 0.7 0.9 0.995];

% white & poldrack use these quantiles
qntls = [0 .2, .4, .6, .8, .95];
%qntls = [0 1/3 2/3 1];

% for each subject and previous response
subjects    = unique(data.subj_idx)';
sessions    = unique(data.session)';

CRF         = nan(numel(sessions), numel(subjects), 2, numel(qntls)-1);
CRF_norm    = nan(numel(sessions), numel(subjects), 1, numel(qntls)-1);

plotedges   = nan(numel(sessions), numel(subjects), numel(qntls)-1);
edgeIx      = nan(numel(sessions), numel(subjects), numel(qntls)-1);

colors      = cbrewer('div', 'PiYG', 7);
colors      = colors([1 end], :);
prevresps   = unique(data.prevresp);

for sj = subjects,
    for s = sessions,
        
        % divide into RT bins without the previous response conditioning
        thisdat = data(data.subj_idx == sj & data.session == s, :);
        
        if ~isempty(thisdat),
            
            % conditional accuracy function!
            edges                       = quantile(thisdat.rt, qntls);
            thisdat.rtbin               = discretize(thisdat.rt, edges);
            thisdat(isnan(thisdat.rtbin), :) = [];
            
            for r = unique(thisdat.rtbin)',
                for p = 1:2,
                    [~, CRF(s, find(sj==subjects), p, r)]        = ...
                      dprime(thisdat.stimulus(thisdat.rtbin == r & thisdat.prevresp == prevresps(p)), ...
                        thisdat.response(thisdat.rtbin == r & thisdat.prevresp == prevresps(p)));
                end
                [~, CRF_norm(s, find(sj==subjects), 1, r)]   = ...
                    dprime(thisdat.stimulus(thisdat.rtbin == r), thisdat.response(thisdat.rtbin == r));
                plotedges(s, find(sj==subjects), r)     = nanmean(thisdat.rt(thisdat.rtbin == r));
                edgeIx(s, find(sj==subjects), r)        = nanmean(thisdat.rtbin(thisdat.rtbin == r));
            end
            assert(all(diff(plotedges(s, find(sj==subjects), :)) > 0), 'plotedges does not increase');
        end
    end
end

% make criterion go in the direction of % choices
CRF = -CRF;
CRF_norm = -CRF_norm;

% normalise
CRF = bsxfun(@minus, CRF, CRF_norm);

% average over sessions
CRF         = nanmean(CRF, 1);
plotedges   = nanmean(plotedges, 1);
edgeIx      = nanmean(edgeIx, 1);

if numel(subjects) > 1,
    h = boundedline(squeeze(nanmean(edgeIx, 2)), squeeze(nanmean(CRF, 2)), ...
        permute(squeeze(nanstd(CRF, [], 2)) ./ sqrt(numel(unique(data.subj_idx))), [2 3 1]), ...
        'cmap', colors, 'alpha');
    
    [s, pval] = ttest(squeeze(CRF(:, :, 1, :)), squeeze(CRF(:, :, 2, :)));
    plot(find(s == 1), -0.3 * ones(size(find(s==1))), '.k', 'markersize', 10);
    
else
    hold on;
    for c = 1:2,
        h = plot(squeeze(nanmean(edgeIx, 2)), ...
            squeeze(nanmean(CRF(:, :, c, :), 2)), 'color', colors(c, :));
    end
end

xlabel('RT quantile'); ylabel('Bias towards choice 1')
axisNotSoTight;
if numel(subjects) > 1,
     ylim([-0.5 0.5]);
end
set(gca, 'xtick', 1:length(nanmean(edgeIx)), 'xticklabel', qntls(2:end));
end

% ========================================== %
% function for fixed-effects RT distributions
% ========================================== %

function h = plotRTdist(data)

colors      = cbrewer('div', 'PiYG', 7);
colors      = colors([1 end], :);

prevresp = [-1 1];
hold on;
for p = 1:2,
    h1{p} =  histogram(data.rt(data.prevresp == prevresp(p)), ...
        'displaystyle', 'stairs', 'edgecolor', colors(p, :), 'normalization', 'pdf', ...
        'binwidth', 0.05, 'linewidth', 1);
end

h = [h1{:}];
end