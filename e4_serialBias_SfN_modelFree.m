function e4_serialBias_SfN_modelFree

addpath(genpath('~/code/Tools'));
warning off; close all; clear;
global datasets datasetnames
datasets = {'RT_RDK', 'projects/0/neurodec/Data/MEG-PL', 'Anke_2afc_sequential'};
datasetnames = {'RT', '2IFC', 'Anke'};

set(groot, 'defaultaxesfontsize', 7, 'defaultaxestitlefontsizemultiplier', 1.1, ...
    'defaultaxestitlefontweight', 'bold', ...
    'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters');

% ========================================== %
% MODELFREE MEASURE OF BIAS
% RT DISTRIBUTIONS
% ========================================== %
cnt = 1;

for d = 1:length(datasets),
    
    % load data
    clearvars -except d datasets cnt;
    csvfile = dir(sprintf('~/Data/%s/HDDM/*.csv', datasets{d}));
    csvfile = csvfile(arrayfun(@(x) ~strcmp(x.name(1),'.'), csvfile)); % remove hidden files
    alldata = readtable(sprintf('~/Data/%s/HDDM/%s', datasets{d}, csvfile.name));
    
    % define repeaters and alternators based on dc
    dat = readtable(sprintf('~/Data/%s/HDDM/summary/allindividualresults.csv', ...
        datasets{d}));
    if d == 3,
        alldata = alldata(alldata.transitionprob == 0.5, :);
        alldata.stimulus = sign(alldata.motionenergy);
        groups{1} = dat.subjnr(dat.v_prevresp_neutral__regressdcprevresp < 0);
        groups{2} = dat.subjnr(dat.v_prevresp_neutral__regressdcprevresp > 0);
    else
        groups{1} = dat.subjnr(dat.v_prevresp__regressdcprevrespstim < 0);
        groups{2} = dat.subjnr(dat.v_prevresp__regressdcprevrespstim > 0);
    end
    
    for gr = 1:length(groups),
        subplot(4,4,cnt); cnt = cnt + 1;
        plotRTdistributions(alldata(ismember(alldata.subj_idx, ...
            groups{gr}) & alldata.stimulus == -1, :));
        % stim == 1
        
        subplot(4,4,cnt); cnt = cnt + 1;
        % stim == -1
        plotRTdistributions(alldata(ismember(alldata.subj_idx, ...
            groups{gr}) & alldata.stimulus == 1, :));
        
    end
end

print(gcf, '-dpdf', '~/Data/serialHDDM/fig2_driftrate.pdf');

end

function plotRTdistributions(data)

colors = cbrewer('div', 'BrBG', 5);
colors = colors([1 end], :);
ylims = [0 3];
data.rt(data.response == 0) = -data.rt(data.response == 0);

hold on;
prevchoices = unique(data.prevresp);
prevchoices(isnan(prevchoices)) = [];

edges       = -3:0.05:3;
% normalisation of the responses?
distFun     = @(x) histcounts(x, edges, 'normalization', 'pdf');
plotedges   = mean([edges', circshift(edges, [0 1])'], 2);
plotedges   = plotedges(2:end);

%%%%%%%%%%%%%%
for p = 1:length(prevchoices),
    RTs     = nan(max(data.subj_idx), numel(plotedges));
    medians = nan(max(data.subj_idx), 2);
    
    for sj = unique(data.subj_idx)',
        RTs(sj, :)      = distFun(data.rt(data.subj_idx == sj & data.prevresp == prevchoices(p)));
        medians(sj, 1)  = nanmedian(data.rt(data.subj_idx == sj & ...
            data.prevresp == prevchoices(p) & data.rt < 0));
        medians(sj, 2)  = nanmedian(data.rt(data.subj_idx == sj & ...
            data.prevresp == prevchoices(p) & data.rt > 0));
    end
    
    % for both positive and negative RTs
    plot([nanmean(medians(:, 1)) nanmean(medians(:, 1))], [ylims(2)*0.8 ylims(2)], 'color', colors(p, :), 'linewidth', 0.5);
    plot([nanmean(medians(:, 2)) nanmean(medians(:, 2))], [ylims(2)*0.8 ylims(2)], 'color', colors(p, :), 'linewidth', 0.5);

    % plot
    boundedline(plotedges, squeeze(nanmean(RTs)), permute(squeeze(nanstd(RTs)) ...
        ./ sqrt(numel(unique(data.subj_idx))), [2 3 1]), 'cmap', colors(p, :), 'alpha');
end

xlabel('Response time (s)'); ylabel('Probability');
xlim([-3 3]); box off;
ylim(ylims);

end