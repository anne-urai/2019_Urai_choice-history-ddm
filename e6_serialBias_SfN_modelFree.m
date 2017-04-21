function e6_serialBias_SfN_modelFree

addpath(genpath('~/code/Tools'));
warning off; close all; clear;
global datasets datasetnames
datasets = {'RT_RDK', 'projects/0/neurodec/Data/MEG-PL', 'NatComm', 'Anke_2afc_neutral'};
datasetnames = {'RT', '2IFC', 'NatComm', 'Anke_neutral'};

set(groot, 'defaultaxesfontsize', 7, 'defaultaxestitlefontsizemultiplier', 1, ...
    'defaultaxestitlefontweight', 'normal', ...
    'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters');

% ========================================== %
% MODELFREE MEASURE OF BIAS
% RT DISTRIBUTIONS
% ========================================== %

for d = 1:length(datasets),
    
    % load data
    clearvars -except d datasets cnt datasetnames whichSJ ylims;
    csvfile = dir(sprintf('~/Data/%s/HDDM/*.csv', datasets{d}));
    csvfile = csvfile(arrayfun(@(x) ~strcmp(x.name(1),'.'), csvfile)); % remove hidden files
    alldata = readtable(sprintf('~/Data/%s/HDDM/%s', datasets{d}, csvfile.name));
    
    if d > 2, % in case of natcomm and anke's datasets
        alldata.stimulus = sign(alldata.motionenergy);
    end
    
    subjects = unique(alldata.subj_idx)';
    for sj = subjects,
        
        clf; cnt = 1;
        stimuli = [-1 1];
        for s = 1:2,
            subplot(4,4,cnt); cnt = cnt + 4;
            h = plotRTdistributions(alldata(alldata.subj_idx == sj ...
                & alldata.stimulus == stimuli(s),:));
            xlim([-3 3]);
            if s == 2, xlabel('Response time (s)');
                title(sprintf('stimulus %d', stimuli(s)));
            else
                title(sprintf('stimulus %d, P%02d', stimuli(s), sj));
            end
        end
        
        % layout and save
        l = legend(h, {'prevresp -1', 'prevresp 1'}, 'Location', 'South');
        l.Box = 'off';
        l.Position(2) = l.Position(2) - 0.12;
        print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/fig4_RTdistributions_%s_P%02d.pdf', datasetnames{d}, sj));
    end
end
end

function h = plotRTdistributions(data)

colors = cbrewer('div', 'PiYG', 7);
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

%     whichFun = 'twosided';
%     switch whichFun
%         case 'kde'
%             distFun     = @(x) ksdensity(x, edges);
%             %[X,fVal,exitFlag,solverOutput] = exgauss_fit(y);
%         case 'egfit'
%             % from http://www.psy.ulaval.ca/?pid=1529
%             distFun     = @(x) egfit(x(~isnan(x)));
%         case 'twosided'
%             distFun     = @(x) ksdensity(x, edges);
%             edges       = [-fliplr(edges) edges];
%     end
%
%     subjects    = unique(data.subj_idx);
%     RTs         = nan(numel(subjects), 2, numel(edges));
%     params      = nan(numel(subjects), 2, 3);
%
%     for p = 1:length(reps),
%
%         for sj = unique(data.subj_idx)',
%             sjidx = find(sj == subjects);
%
%             switch whichFun
%                 case 'egfit'
%                     % separately for positive and negative distribution, compute pdf
%                     a1 = distFun(data.rt(data.repeat == reps(p) & data.subj_idx == sj));
%
%                     % extract the parameters of the ex-gaussian distribution
%                     params(sjidx, p, :)   = a1;
%                     RTs(sjidx, p, :)      = exgauss_pdf(edges, a1);
%                 case 'kde';
%                     RTs(sjidx, p, :)      = distFun(data.rt(data.repeat == reps(p) & data.subj_idx == sj));
%                 case 'twosided'
%                     stimuli = [-1 1];
%                     for s = 1:length(stimuli),
%                         tmp{s} = distFun(data.rt(data.repeat == reps(p) & data.stimulus == stimuli(s) & data.subj_idx == sj));
%                     end
%                     RTs(sjidx, p, :) = cat(2, fliplr(tmp{1}), tmp{2});
%             end
%         end
%     end
%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % compute average over subjects
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     for p = 1:2,
%         % compute the mean across the group
%         switch whichFun
%             case 'egfit'
%                 RTs_group = exgauss_pdf(edges, squeeze(mean(params(:, p, :), 1)));
%             otherwise
%                 RTs_group = squeeze(nanmean(RTs(:, p, :), 1));
%         end
%
%         h(p) = plot(edges, RTs_group, 'color', colors(p, :));
%     end
%
%     xlabel('Response time (s)'); ylabel('Fraction of trials');
%     axis tight; box off;

end