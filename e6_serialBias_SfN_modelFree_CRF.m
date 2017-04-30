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
if 0,
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
            h = plotRTdistributions(alldata(ismember(alldata.subj_idx, groups{gr}), :));
            title({datasetnames{d}, sprintf('%s n = %d', groupnames{gr}, numel(groups{gr}))});
            axis square;
        end
        
    end
    
    %suplabel(['Data: ' datasetnames{d}], 't', [sph{:}]);
    l               = legend(h, {'Prev choice -1', 'Prev choice 1'}, 'Location', 'South');
    l.Box           = 'off';
    l.Position(2)   = l.Position(2) - 0.15;
    
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/fig4_CRF.pdf'));
end

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
        h = plotRTdistributions(alldata(alldata.subj_idx == sj, :));
        title(sprintf('P%02d', sj));
        axis square;
    end
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/fig4_CRF_data%d_allsubjects.pdf', d));
end

end

function h = plotRTdistributions(data)

qntls = [0.1 0.3 0.5 0.7 0.9 0.995];

% for each subject and previous response
subjects    = unique(data.subj_idx)';
sessions    = unique(data.session)';

CRF         = nan(numel(sessions), numel(subjects), 2, numel(qntls)-1);
CRF_norm    = nan(numel(sessions), numel(subjects), 1, numel(qntls)-1);

plotedges   = nan(numel(sessions), numel(subjects), numel(qntls)-1);
edgeIx      = nan(numel(sessions), numel(subjects), numel(qntls)-1);

colors      = cbrewer('div', 'PiYG', 7);
colors      = colors([1 end], :);

for sj = subjects,
    for s = sessions,
        
        % divide into RT bins without the previous response conditioning
        thisdat = data(data.subj_idx == sj & data.session == s, :);
        
        if ~isempty(thisdat),
            % conditional accuracy function!
            edges                       = quantile(thisdat.rt, qntls);
            thisdat.rtbin               = discretize(thisdat.rt, edges);
            [gr, prevrespidx, rtidx]    = findgroups(thisdat.prevresp, thisdat.rtbin);
            
            tmp_c                                 = splitapply(@nanmean, thisdat.response, gr);
            CRF(s, find(sj==subjects), 1, :)      = tmp_c(1:numel(qntls)-1);
            CRF(s, find(sj==subjects), 2, :)      = tmp_c(numel(qntls):end);
            
            % check how well the division worked
            plotedges(s, find(sj==subjects), :)   = splitapply(@nanmean, thisdat.rt, findgroups(thisdat.rtbin));
            edgeIx(s, find(sj==subjects), :)      = splitapply(@nanmean, thisdat.rtbin, findgroups(thisdat.rtbin));
            CRF_norm(s, find(sj==subjects), 1, :) = splitapply(@nanmean, thisdat.response, findgroups(thisdat.rtbin));
        end
    end
end

% normalise
CRF = bsxfun(@minus, CRF, CRF_norm);

% average over sessions
CRF         = nanmean(CRF, 1);
plotedges   = nanmean(plotedges, 1);
edgeIx      = nanmean(edgeIx, 1);

if numel(subjects) > 1,
    h = boundedline(squeeze(nanmean(plotedges, 2)), squeeze(nanmean(CRF, 2)), permute(squeeze(nanstd(CRF, [], 2)) ...
        ./ sqrt(numel(unique(data.subj_idx))), [3 1 2]), 'cmap', colors, 'alpha');
else
    hold on;
    for c = 1:2,
        h = plot(squeeze(nanmean(plotedges, 2)), ...
            squeeze(nanmean(CRF(:, :, c, :), 2)), 'color', colors(c, :));
    end
end
% hold on;
% for c = 1:2,
%     h = ploterr(nanmean(plotedges), squeeze(nanmean(CRF(:, c, :))), ...
%         nanstd(plotedges) ./ sqrt(length(subjects)), ...
%         squeeze(nanstd(CRF(:, c, :))) ./ sqrt(length(subjects)), 'abshhxy', 0);
%     set(h(1), 'color', colors(c, :));
%     set(h(2), 'color', colors(c, :));
%     set(h(3), 'color', colors(c, :));
% end

xlabel('RT (ms)'); ylabel('P(choice = 1)')
axisNotSoTight;
if numel(subjects) > 1,
    ylim([-0.10 0.10]);
end
try
    set(gca, 'xtick', squeeze(nanmean(plotedges, 2)));
end

end