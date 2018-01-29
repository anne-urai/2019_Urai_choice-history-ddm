function e6_serialBias_SfN_modelFree

addpath(genpath('~/code/Tools'));
warning off; close all; clear;
global datasets datasetnames
datasets = {'RT_RDK', 'projects/0/neurodec/Data/MEG-PL', 'NatComm', 'Anke_2afc_neutral'};
datasetnames = {'RT', '2IFC', 'NatComm', 'Anke neutral'};

set(groot, 'defaultaxesfontsize', 8, 'defaultaxestitlefontsizemultiplier', 1, ...
    'defaultaxestitlefontweight', 'normal', ...
    'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters');

% ========================================== %
% MODELFREE MEASURE OF BIAS
% RT DISTRIBUTIONS
% ========================================== %

for d = 1:length(datasets),
    clearvars -except d datasets cnt datasetnames whichSJ ylims;
    
    % load individual serial choice bias
    dat = readtable(sprintf('~/Data/%s/HDDM/summary/allindividualresults.csv', ...
        datasets{d}));
    serialbias = dat.v_prevresp__regressdcprevrespstim(dat.session == 0);
    
    % load data
    csvfile = dir(sprintf('~/Data/%s/HDDM/*.csv', datasets{d}));
    csvfile = csvfile(arrayfun(@(x) ~strcmp(x.name(1),'.'), csvfile)); % remove hidden files
    alldata = readtable(sprintf('~/Data/%s/HDDM/%s', datasets{d}, csvfile.name));
    
    if d > 2, % in case of natcomm and anke's datasets
        alldata.stimulus = sign(alldata.motionenergy);
    end
    
    subjects = unique(alldata.subj_idx)';
    [~, idx] = sort(serialbias); idx = idx';
    for sj = idx,
        
        clf; cnt = 1;
        stimuli = [-1 1];
        for s = 1:2,
            subplot(4,4,cnt); cnt = cnt + 4;
            h = plotRTdistributions(alldata(alldata.subj_idx == subjects(sj) ...
                & alldata.stimulus == stimuli(s),:));
            xlim([-3 3]);
            if s == 2, xlabel('Response time (s)');
                title(sprintf('stimulus %d', stimuli(s)));
            else
                title(sprintf('stimulus %d', stimuli(s)));
                suplabel(sprintf('%s, P%02d, bias %.3f', datasetnames{d}, subjects(sj), serialbias(sj)), 't');
            end
        end
        
        % layout and save
        l = legend(h, {'prevresp -1', 'prevresp 1'}, 'Location', 'South');
        l.Box = 'off';
        l.Position(2) = l.Position(2) - 0.12;
        print(gcf, '-depsc', sprintf('~/Data/serialHDDM/fig4_RTdistributions_%s_P%02d.eps', ...
            datasetnames{d}, sj));
    end
end
end

function h = plotRTdistributions(data)

colors  = cbrewer('div', 'PiYG', 7);
colors1 = colors([1 end], :);
colors2 = colors([2 end-1], :);
hold on;
% step size
edges       = 0:0.03:3;
prevresps   = [-1 1];
resps       = [0 1];

% two functions
% distFun     = @(x) ksdensity(x, edges);
distFun     = @(x) histcounts(x, [edges 3.001], 'normalization', 'pdf');

fitFun      = @(x) exgauss_pdf(edges, egfit(x(~isnan(x))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% split by previous response
% show both distribution and KDE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% do separate fits for choice == 1 and choice == -1 (error and correct,
% depending on the stimulus)
for p = 1:length(prevresps),
    
    for c = 1:length(resps),
        thisdat = data(data.prevresp == prevresps(p) & data.response == resps(c), :);
        dist{p, c}    = distFun(thisdat.rt);
        fitted{p, c}  = fitFun(thisdat.rt);
    end
end

% now first plot the RT distributions by themselves
for p = 1:2,
    for c = 1:2,
        if resps(c) == 0,
            plot(fliplr(-edges), fliplr(dist{p,c}), 'color', colors2(p, :));
        else
            plot(edges, dist{p,c}, 'color', colors2(p, :));
        end
    end
    
end

% now the thicker fits on top
for p = 1:2,
    for c = 1:2,
        if resps(c) == 0,
            plot(fliplr(-edges), fliplr(fitted{p,c}), 'color', colors1(p, :), 'linewidth', 1);
        else
            h(p) = plot(edges, fitted{p,c}, 'color', colors1(p, :), 'linewidth', 1);
        end
    end
end

ylabel('Fraction of trials');
vline(0);
text(-2.5, mean(get(gca, 'ylim'))*1.5, 'choice -1', 'fontsize', 5);
text(1.5, mean(get(gca, 'ylim'))*1.5, 'choice 1', 'fontsize', 5);

end