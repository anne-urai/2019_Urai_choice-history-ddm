function e6_serialBias_SfN_modelFree

addpath(genpath('~/code/Tools'));
warning off; close all; clear;
global datasets datasetnames
datasets = {'RT_RDK', 'projects/0/neurodec/Data/MEG-PL', 'Anke_2afc_neutral', 'NatComm'};
datasetnames = {'RT', '2IFC', 'Anke neutral', 'NatComm'};
ylims = {[0 0.002], [0 0.003], [0 0.004], [0 0.004]};

set(groot, 'defaultaxesfontsize', 7, 'defaultaxestitlefontsizemultiplier', 1.1, ...
    'defaultaxestitlefontweight', 'bold', ...
    'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters');
whichSJ = 'subgroups';

% ========================================== %
% MODELFREE MEASURE OF BIAS
% RT DISTRIBUTIONS
% ========================================== %
clf;
cnt = 1;
for d = 1:length(datasets),
    
    % load data
    clearvars -except d datasets cnt datasetnames whichSJ ylims;
    csvfile = dir(sprintf('~/Data/%s/HDDM/*.csv', datasets{d}));
    csvfile = csvfile(arrayfun(@(x) ~strcmp(x.name(1),'.'), csvfile)); % remove hidden files
    alldata = readtable(sprintf('~/Data/%s/HDDM/%s', datasets{d}, csvfile.name));
    % add repeated vs alternated response
    alldata.repeat = [NaN; (diff(alldata.response) == 0)];
    
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
    
    groupnames = {'alternators', 'repeaters'};
    % call the plotting function
    for gr = 1:length(groups),
        subplot(4,4,cnt); cnt = cnt + 1;
        if numel(groups{gr}) > 0,
            h = plotRTdistributions(alldata(ismember(alldata.subj_idx, ...
                groups{gr}), :), [0 0.005]);
            title(sprintf('%s, %s', groupnames{gr}, datasetnames{d}));
        else
            axis off;
        end
    end
    
end

l = legend(h, {'switch', 'repeat'}, 'Location', 'South');
l.Box = 'off';
l.Position(2) = l.Position(2) - 0.15;
l.Position(1) = l.Position(1) - 0.1;

print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/fig4_RTdistributions_%d.pdf', d));
end

function h = plotRTdistributions(data, ylims)

colors = cbrewer('div', 'PuOr', 7);
colors = colors([2 end-1], :);
hold on;
reps = [0 1];

% step size
edges       = 0:0.01:3;

whichFun = 'twosided';
switch whichFun
    case 'kde'
        distFun     = @(x) ksdensity(x, edges);
        %[X,fVal,exitFlag,solverOutput] = exgauss_fit(y);
    case 'egfit'
        % from http://www.psy.ulaval.ca/?pid=1529
        distFun     = @(x) egfit(x(~isnan(x)));
    case 'twosided'
        distFun     = @(x) ksdensity(x, edges);
        edges       = [-fliplr(edges) edges];
end

%%%%%%%%%%%%%%
subjects    = unique(data.subj_idx);
RTs         = nan(numel(subjects), 2, numel(edges));
params      = nan(numel(subjects), 2, 3);

for p = 1:length(reps),
    
    for sj = unique(data.subj_idx)',
        sjidx = find(sj == subjects);
        
        switch whichFun
            case 'egfit'
                % separately for positive and negative distribution, compute pdf
                a1 = distFun(data.rt(data.repeat == reps(p) & data.subj_idx == sj));
                
                % extract the parameters of the ex-gaussian distribution
                params(sjidx, p, :)   = a1;
                RTs(sjidx, p, :)      = exgauss_pdf(edges, a1);
            case 'kde';
                RTs(sjidx, p, :)      = distFun(data.rt(data.repeat == reps(p) & data.subj_idx == sj));
            case 'twosided'
                stimuli = [-1 1];
                for s = 1:length(stimuli),
                    tmp{s} = distFun(data.rt(data.repeat == reps(p) & data.stimulus == stimuli(s) & data.subj_idx == sj));
                end
                RTs(sjidx, p, :) = cat(2, fliplr(tmp{1}), tmp{2});
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute average over subjects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for p = 1:2,
    % compute the mean across the group
    switch whichFun
        case 'egfit'
            RTs_group = exgauss_pdf(edges, squeeze(mean(params(:, p, :), 1)));
        otherwise
            RTs_group = squeeze(nanmean(RTs(:, p, :), 1));
    end
    
    h(p) = plot(edges, RTs_group, 'color', colors(p, :));
end

xlabel('Response time (s)'); ylabel('Fraction of trials');
axis tight; box off;

end