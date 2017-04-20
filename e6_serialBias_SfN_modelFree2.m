function e6_serialBias_SfN_modelFree2

addpath(genpath('~/code/Tools'));
warning off; close all; clear;
global datasets datasetnames
datasets = {'RT_RDK', 'projects/0/neurodec/Data/MEG-PL', 'Anke_2afc_neutral', 'NatComm'};
datasetnames = {'RT', '2IFC', 'Anke neutral', 'NatComm'};

set(groot, 'defaultaxesfontsize', 7, 'defaultaxestitlefontsizemultiplier', 1.1, ...
    'defaultaxestitlefontweight', 'bold', ...
    'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters');

whichSJ = 'subgroups';

% ========================================== %
% MODELFREE MEASURE OF BIAS
% RT DISTRIBUTIONS
% ========================================== %
clf; cnt = 1;

for d = 1:length(datasets),
    
    % load data
    clearvars -except d datasets cnt datasetnames whichSJ ylims;
    csvfile = dir(sprintf('~/Data/%s/HDDM/*.csv', datasets{d}));
    csvfile = csvfile(arrayfun(@(x) ~strcmp(x.name(1),'.'), csvfile)); % remove hidden files
    alldata = readtable(sprintf('~/Data/%s/HDDM/%s', datasets{d}, csvfile.name));
    
    if d == 4,
        % HDDM still running, get fruend weights
        load('/Users/anne/Data/pupilUncertainty_FigShare/Data/GrandAverage/historyweights_plain.mat');
        
        switch whichSJ
            case 'extremes'
                [~, idx]    = min(dat.response(:, 1));
                groups{1}   = idx;
                [~, idx]    = max(dat.response(:, 1));
                groups{2}   = idx;
            case 'subgroups'
                groups{1}   = find(dat.response(:, 1) < 0);
                groups{2}   = find(dat.response(:, 1) > 0);
        end
        alldata.stimulus = sign(alldata.motionenergy);
    else
        
        % define repeaters and alternators based on dc
        dat = readtable(sprintf('~/Data/%s/HDDM/summary/allindividualresults.csv', ...
            datasets{d}));
        if d == 3,
            alldata.stimulus = sign(alldata.motionenergy);
            alldata.rt = alldata.rt + 0.5;
        end
        
        switch whichSJ
            case 'extremes'
                % find alternators and repeaters
                [~, idx]    = min(dat.v_prevresp__regressdcprevrespstim);
                groups{1}   = dat.subjnr(idx);
                [~, idx]    = max(dat.v_prevresp__regressdcprevrespstim);
                groups{2}   = dat.subjnr(idx);
            case 'subgroups'
                % find alternators and repeaters
                groups{1}   = find(dat.v_prevresp__regressdcprevrespstim < 0);
                groups{2}   = find(dat.v_prevresp__regressdcprevrespstim > 0);
        end
    end
    
    groupnames  = {'alternators', 'repeaters'};
    for gr = 1:length(groups),
        sph{gr} = subplot(4,4,cnt); cnt = cnt + 1;
        try
            h = plotRTdistributions(alldata(ismember(alldata.subj_idx, ...
                groups{gr}), :));
            title(sprintf('%s, %s', groupnames{gr}, datasetnames{d}));
            axis square;
        end
    end
    
    
end
%suplabel(['Data: ' datasetnames{d}], 't', [sph{:}]);
l = legend(h, {'Prev choice -1', 'Prev choice 1'}, 'Location', 'South');
l.Box = 'off';
l.Position(2) = l.Position(2) - 0.15;
%l.Position(1) = l.Position(1) - 0.1;

print(gcf, '-depsc', sprintf('~/Data/serialHDDM/fig4_CRF.eps'));

end

function h = plotRTdistributions(data)
% conditional accuracy function!
nbins = 4;
edges = quantile(data.rt, nbins-1);

hold on;
binIdx      = @(x) discretize(x, [0 edges 10]);
plotedges   = splitapply(@nanmean, data.rt, binIdx(data.rt));
prevchoices = [-1 1];

CRF         = nan(max(data.subj_idx), 2, nbins);
fullCRF     = nan(max(data.subj_idx), 1, nbins);

colors      = cbrewer('div', 'PRGn', 7);
colors      = colors([1 end], :);

binIdx      = @(x) discretize(x, [0 quantile(x, nbins-1) 10]);

% for each subject and previous response
for sj = unique(data.subj_idx)',
    
    % divide into RT bins without the previous response conditioning
    thisdat = data(data.subj_idx == sj, :);
    thisdat.binIdx  = binIdx(thisdat.rt);
    
    for p = 1:length(prevchoices),
        % split by rt
        for e = 1:nbins,
            CRF(sj, p, e) = nanmean(thisdat.response(thisdat.binIdx == e & thisdat.prevresp == prevchoices(p))); ...
        end
    end
    
    % also compute one on regardless of prevchoice
    for e = 1:nbins,
        fullCRF(sj, 1, e) = nanmean(thisdat.response(thisdat.binIdx == e)); ...
    end

end
plot(plotedges([1 end]),[0 0],  'color', [0.5 0.5 0.5]);

% normalize
CRF = bsxfun(@minus, CRF, fullCRF);

% plot
h = boundedline(plotedges, squeeze(nanmean(CRF)), permute(nanstd(CRF) ...
    ./ sqrt(numel(unique(data.subj_idx))), [3 1 2]), 'cmap', colors, 'alpha');

xlabel('RT quantiles'); ylabel('P(choice == 1)')
% ylim([-0.1 0.1]); set(gca, 'ytick', [-0.1 0 0.1]);
xlim([min(plotedges)*0.8 max(plotedges)]);
axis square;

end